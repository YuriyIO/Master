git clone https://gitlab.inria.fr/scotch/scotch.git

cd scotch/ &&  mkdir build && cd build/

cmake .. -DCMAKE_INSTALL_PREFIX=/home_edu/edu-cmc-sqi23/edu-cmc-sqi23-25/local_libs/scotch     -DCMAKE_C_FLAGS="-D_SCOTCHyylex=yylex -D_SCOTCHyyparse=yyparse -D_SCOTCHyyerror=yyerror -D_SCOTCHyylval=yylval -D_SCOTCHyychar=yychar -D_SCOTCHyydebug=yydebug -D_SCOTCHyynerrs=yynerrs"

make -j1

make install

ln -s lib64 lib

cd ..
