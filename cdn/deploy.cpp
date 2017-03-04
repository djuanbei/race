#include "deploy.h"

#include "lib_io.h"
#include "lib_time.h"
#include<vector>
#include<stdlib.h>
#include <stdio.h>
#include<string.h>


using namespace std;
//You need to complete the function 
void deploy_server(char * topo[MAX_EDGE_NUM], int line_num,char * filename)
{
  
  char line[10000  ];
  strcpy(line, topo[ 0 ]);
  char *temp=  strtok( line, " " );
  int V=atoi(temp  );
  temp=  strtok( NULL, " " );
  int E=atoi(temp  );
  temp=  strtok( NULL, " " );
  int D= atoi( temp );

  strcpy(line, topo[ 2 ]);
  temp=  strtok( line, " " );
  
  int serice_value=atoi( temp );
  
  
  int i=4;
  strcpy(line, topo[ 4 ]);
  vector<int> srcs, snks, caps, ws;
  while( strlen( line )>0 ){
    int src, snk, cap, w;
    
    temp=  strtok( line, " " );
    src=atoi(temp  );
    temp=  strtok( NULL, " " );
    snk=atoi(temp  );
    temp=  strtok( NULL, " " );
    cap=atoi(temp  );
    temp=  strtok( NULL, " " );
    w=atoi(temp  );

    srcs.push_back( src );
    snks.push_back( snk );
    caps.push_back( cap );
    ws.push_back( w );
    
    strcpy(line, topo[ i++ ]);
  }
  
  directed_graph graph;
  graph.initial(srcs, snks, ws);
  
  
	// Output demo``''
  char * topo_file = (char *)"17\n\n0 8 0 20\n21 8 0 20\n9 11 1 13\n21 22 2 20\n23 22 2 8\n1 3 3 11\n24 3 3 17\n27 3 3 26\n24 3 3 10\n18 17 4 11\n1 19 5 26\n1 16 6 15\n15 13 7 13\n4 5 8 18\n2 25 9 15\n0 7 10 10\n23 24 11 23";

  write_result(topo_file, filename);

}



