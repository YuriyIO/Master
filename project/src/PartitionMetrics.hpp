#ifndef PARTITION_METRICS_HPP
#define PARTITION_METRICS_HPP

#include "CSR.hpp"
#include "Partition.hpp"
#include <vector>
#include <string>

class PartitionMetrics {
public:
    const std::string outputFolder;
    const std::string vizFolder;

    std::vector<int> ghost_global_ids;
    std::vector<int> ghost_part_ids;

    std::vector<int> verticesCounter;
    std::vector<double> vertexImbalance;
    std::vector<long long> edgesCounter;
    std::vector<double> edgeImbalance;
    std::vector<int> boundaryPerPart;
    std::vector<int> componentsPerPart;

    std::vector<std::vector<int>> neighbours;

    PartitionMetrics(const std::string outputFolder, const std::string vizFolder);

    void save_partition_info(const Partition& partition, const int actualParts, const CSR& csr);

    void save_graph_topology(const CSR& csr);
    void save_partition_vec(const Partition& partition, const CSR& csr);

private:
    struct Edge {
        int p1, p2;
        bool operator<(const Edge& other) const;
        bool operator==(const Edge& other) const;
    };

    bool graph_structure_initialized = false;
    std::vector<int> send_counts, recv_counts;
    std::vector<int> send_displs, recv_displs; 
    std::vector<int> export_indices;
    std::vector<int> sorting_p;
    int total_send, total_recv;


    int findOwnerRank(const int v, const std::vector<int>& vtxdist);
    int getPartitionOfVertex(int v, const CSR& csr, const std::vector<int>& partition);
    void exchange_ghost_parts(const std::vector<int>& partition, const int actualParts, const CSR& csr);
    void exchange_ghost_parts2(const std::vector<int>& partition, const int actualParts, const CSR& csr);
    void computeVertexImbalance(const std::vector<int>& partition, const int actualParts);
    void computeEdgeImbalance(const std::vector<int>& partition, const int actualParts, const CSR& csr);
    void computeBoundaryVertices(const std::vector<int>& partition, const int actualParts, const CSR& csr);
    void countCutEdges(const std::vector<int>& partition, const int actualParts, std::vector<int>* cutEdgesPerPart, int& totalCutEdges, const CSR& csr);
    void countNeighbours(const std::vector<int>& partition, const int actualParts, const CSR& csr);
    std::string getFileName(const std::string& filename);
    void saveMetricsToFile(const Partition& partition, const int cutEdge, const int rank, const std::string filename);
    void computeConnectedComponents(const std::vector<int>& partition, const int actualParts, const CSR& csr);
};

#endif