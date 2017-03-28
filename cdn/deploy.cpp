#include "deploy.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include "lib_io.h"
#include "lib_time.h"


#include "locchoose.h"

using namespace std;

using namespace raptor;


// You need to complete the function
void deploy_server(char *topo[MAX_EDGE_NUM], int line_num, char *filename) {
  int max_line_length = 1024;
  size_t str_len = 1023;
  char *delim = " ";

  char *temp_str = strtok(topo[0], delim);
  int network_node_num = atoi(temp_str);
  temp_str = strtok(NULL, delim);
  int network_link_num = atoi(temp_str);
  temp_str = strtok(NULL, delim);
  int user_node_num = atoi(temp_str);

  int vir_source = network_node_num;
  int vir_target = vir_source + 1;

  vector<int> network_node_user_map(network_node_num, -100);

  temp_str = strtok(topo[2], delim);

  int server_price = atoi(temp_str);

  vector<int> srcs, snks, caps, weights;

  for (int v = 0; v < network_node_num; v++) {
    srcs.push_back(vir_source);
    snks.push_back(v);
    caps.push_back(0);
    weights.push_back(0);
  }

  int line_id = 4;

  int i = 0;
  while (i < network_link_num) {
    i++;
    int src, snk, cap, w;

    temp_str = strtok(topo[line_id++], delim);
    src = atoi(temp_str);
    temp_str = strtok(NULL, delim);
    snk = atoi(temp_str);
    temp_str = strtok(NULL, delim);
    cap = atoi(temp_str);
    temp_str = strtok(NULL, delim);
    w = atoi(temp_str);

    srcs.push_back(src);
    snks.push_back(snk);
    caps.push_back(cap);
    weights.push_back(w);
  }

  line_id++;

  i = 0;
  int target_link_start = 2 * srcs.size();

  vector<int> user_demand;
  while (i < user_node_num) {
    i++;
    int user_node, network_node, bw;
    temp_str = strtok(topo[line_id++], delim);
    user_node = atoi(temp_str);

    temp_str = strtok(NULL, delim);
    network_node = atoi(temp_str);
    temp_str = strtok(NULL, delim);
    bw = atoi(temp_str);

    network_node_user_map[network_node] = user_node;

    srcs.push_back(network_node);
    snks.push_back(vir_target);
    weights.push_back(0);
    caps.push_back(bw);
    user_demand.push_back(bw);
  }

  raptor::undirected_graph graph;
  graph.initial(srcs, snks, weights);

  raptor::Loc_choose loc_choose(
      graph, network_node_num, user_node_num, server_price, vir_source,
      vir_target, target_link_start, network_node_user_map, user_demand, caps);

  char *topo_file = loc_choose.solve();

  write_result(topo_file, filename);
  delete[] topo_file;
}
