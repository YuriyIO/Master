#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <filesystem>
#include <iomanip>
#include <algorithm>
#include <sstream>

namespace fs = std::filesystem;

struct DetailStats {
    std::string avg = "N/A", diff = "N/A", avg_dev = "N/A", std_dev = "N/A";
    
    std::string toString() const {
        if (avg == "N/A") return "N/A";
        return "avg: " + avg + "  diff: " + diff + "  avg_dev: " + avg_dev + "  std_dev: " + std_dev;
    }
};

struct PartitionReport {
    std::string reportTitle;
    int libraryRank;
    std::string edgeCut = "N/A";
    std::string timeTotal = "N/A";
    std::string timeSolve = "N/A";
    
    std::map<std::string, DetailStats> details;
};

int getLibraryPriority(const std::string& title) {
    if (title.find("ParHIP") != std::string::npos)    return 1;
    if (title.find("dKaMinPar") != std::string::npos) return 2;
    if (title.find("ParMETIS") != std::string::npos)  return 3;
    if (title.find("PT-SCOTCH") != std::string::npos) return 4;
    if (title.find("Zoltan2") != std::string::npos)   return 6;
    if (title.find("Zoltan") != std::string::npos)    return 5;
    return 99;
}

std::string trim(const std::string& s) {
    size_t first = s.find_first_not_of(" = \t\r\n");
    if (first == std::string::npos) return "";
    size_t last = s.find_last_not_of(" = \t\r\n");
    return s.substr(first, (last - first + 1));
}

std::string extractValue(const std::string& line, const std::string& key) {
    size_t pos = line.find(key);
    if (pos == std::string::npos) return "N/A";
    std::string sub = line.substr(pos + key.length());
    std::stringstream ss(sub);
    std::string val;
    ss >> val;
    return val;
}

void parseLineToDetail(const std::string& line, DetailStats& d) {
    d.avg = extractValue(line, "avg:");
    d.diff = extractValue(line, "diff:");
    d.avg_dev = extractValue(line, "avg_dev:");
    d.std_dev = extractValue(line, "std_dev:");
}

int main(int argc, char* argv[]) {
    if (argc < 2) return 1;

    std::string folderPath = argv[1];
    std::vector<PartitionReport> reports;
    std::string fileName;

    for (const auto& entry : fs::directory_iterator(folderPath)) {
        if (entry.is_regular_file() && entry.path().filename().string().find("_info") != std::string::npos) {
            std::ifstream file(entry.path());
            if (!file.is_open()) continue;

            PartitionReport rep;
            std::string line;
            while (std::getline(file, line)) {
                if (line.find("PARTITION REPORT:") != std::string::npos) {
                    rep.reportTitle = trim(line.substr(line.find(":") + 1));
                }
                else if (line.find("File:") != std::string::npos) {
                    std::string value = trim(line.substr(line.find(":") + 1));
                    fileName = std::filesystem::path(value).stem().string();
                }
                else if (line.find("Edge Cut:") != std::string::npos) rep.edgeCut = extractValue(line, "Edge Cut:");
                else if (line.find("Time Total (sec):") != std::string::npos) rep.timeTotal = extractValue(line, "Time Total (sec):");
                else if (line.find("Time Solve Only (sec):") != std::string::npos) rep.timeSolve = extractValue(line, "Time Solve Only (sec):");
                else if (line.find("Vertices:") != std::string::npos) parseLineToDetail(line, rep.details["Vertices"]);
                else if (line.find("Edges:") != std::string::npos) parseLineToDetail(line, rep.details["Edges"]);
                else if (line.find("Boundary:") != std::string::npos) parseLineToDetail(line, rep.details["Boundary"]);
                else if (line.find("Neighbors:") != std::string::npos) parseLineToDetail(line, rep.details["Neighbors"]);
                else if (line.find("Comps:") != std::string::npos) parseLineToDetail(line, rep.details["Comps"]);
            }
            rep.libraryRank = getLibraryPriority(rep.reportTitle);
            reports.push_back(rep);
        }
    }

    std::sort(reports.begin(), reports.end(), [](const PartitionReport& a, const PartitionReport& b) {
        if (a.libraryRank != b.libraryRank) return a.libraryRank < b.libraryRank;
        return a.reportTitle < b.reportTitle;
    });
 
    fs::path outPath = fs::path(folderPath) / (fileName + "_stats.txt");
    std::ofstream out(outPath);
    if (!out.is_open()) return 1;

    std::vector<std::pair<std::string, std::string>> allMetrics = {
        {"Edge Cut", "edge_cut"},
        {"Boundary", "detail"},
        {"Neighbors", "detail"},
        {"Comps", "detail"},
        {"Vertices", "detail"},
        {"Edges", "detail"},
        {"Time: Total, Solve Only (sec)", "time_combined"}
    };

    for (const auto& m : allMetrics) {
        out << "==================== " << m.first << " ====================\n";

        for (size_t i = 0; i < reports.size(); ++i) {
            if (i > 0 && reports[i].libraryRank != reports[i-1].libraryRank) {
                out << "\n";
            }

            out << std::left << std::setw(50) << reports[i].reportTitle << ": ";

            if (m.second == "edge_cut") {
                out << reports[i].edgeCut;
            } else if (m.second == "time_combined") {
                out << reports[i].timeTotal << ", " << reports[i].timeSolve;
            } else {
                out << reports[i].details[m.first].toString();
            }
            out << "\n";
        }
        out << "\n\n";
    }

    out.close();
    return 0;
}