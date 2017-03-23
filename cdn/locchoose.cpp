
#include<string.h>

#include<algorithm>

#include <sys/resource.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>

#include "locchoose.h"


#ifndef CLOCK_MONOTONIC
#define CLOCK_MONOTONIC 1
#endif

namespace raptor {
using namespace std;

static inline double systemTime(void) {
  struct timespec start;
#ifndef __MACH__
  clock_gettime(CLOCK_MONOTONIC, &start);
#endif
  return start.tv_sec + start.tv_nsec / 1000000000.0;
}


void Loc_choose::initial() {
  vector<pair<int, int>> demandC;
  int tempC = totCap;

  for (int network_node = 0; network_node < network_node_num; network_node++) {
    if (network_node_user_map[network_node] > -1) {
      int user_node = network_node_user_map[network_node];

      user_to_network_map[user_node] = network_node;

      int outDegree = graph.getOutDegree(network_node);
      int sumBw = 0;
      for (int i = 0; i < outDegree; i++) {
        int link = graph.getAdj(network_node, i);
        sumBw += orignal_caps[link / 2];
      }
      /**
       * if for a directly connect user node when sum of  out link
       *bandwidth less than user demand then this node must be a server
       *
       */

      if (sumBw < 2 * user_demand[user_node] ) {
        
        user_direct_server[user_node] = true;
        choosedServer.insert(network_node);
        tempC -= user_demand[user_node];
      } else {
        demandC.push_back(make_pair(user_demand[user_node], user_node));
      }
    }
  }

  sort(demandC.rbegin(), demandC.rend());
  for(size_t i=0; i+1< demandC.size(); i++){
    /**
     *  when one user demand greater or equal then sum of all
     * other users demand then this user direct connect node must be
     * a server
     *
     */
    int user_node =demandC[i].second;// it->second;
    int network_node = user_to_network_map[user_node];

    if (!user_direct_server[user_node]) {
      if (2 * (demandC[i].first) >= tempC) {
        user_direct_server[user_node] = true;
        choosedServer.insert(network_node);
        tempC -= demandC[i].first;
      }
    }
  }

  smallest_user = INF;
  large_user = 0;
  for (int user_node = 0; user_node < user_node_num; user_node++) {
    if (user_demand[user_node] < smallest_user) {
      smallest_user = user_demand[user_node];
    }
    if (user_demand[user_node] > large_user) {
      large_user = user_demand[user_node];
    }
  }
}


void Loc_choose::lower_update() {
  

  int cap = 0;
  int temp_value = 0;
  for (size_t i = 0; i < server_candiate.size(); i++) {
    if (cap + server_candiate[i].success_bw < totCap) {
      cap += server_candiate[i].success_bw;
      temp_value += server_candiate[i].total_price;
    } else {
      temp_value += server_candiate[i].total_price *
                    ((totCap - cap) / (server_candiate[i].success_bw + 0.01));
      if (value_lower > temp_value) {
        value_lower = temp_value;
        break;
      }
    }
  }


  vector<vector<int>> cover_domain(network_node_num);

  for (int user_node = 0; user_node < user_node_num; user_node++) {
    for (vector<int>::iterator it = server_candiate_locs[user_node].begin();
         it != server_candiate_locs[user_node].end(); it++) {
      cover_domain[*it].push_back(user_node);
    }
  }

  Fixed_heap Q( network_node_num);

  vector<int> connect_num(network_node_num);

  for (int network_node = 0; network_node < network_node_num; network_node++) {
    int len = cover_domain[network_node].size();
    Q.push(make_pair(-len, network_node));
    connect_num[network_node] = -len;
  }

  
  vector<bool> pass(user_node_num, false);
  int num = user_node_num;
  
  int minServerNum =  INF;

  
  for (int i = 0; i< 10; i++) {
    
    num = user_node_num;
    fill(pass.begin(  ), pass.end(  ), false);
    
    Q.clear(  );
    for (int network_node = 0; network_node < network_node_num;
         network_node++) {
      int len = cover_domain[network_node].size() +raoDong(  );

      Q.push(make_pair(-len, network_node));
      connect_num[network_node] = -len;
    }

    Server_loc tempS1;

    vector<bool> pass(user_node_num, false);
    int num = user_node_num;

    while (!Q.empty()) {
      pair<int,int> p = Q.top();
      int network_node = p.second;
      tempS1.locs.insert(network_node);

      Q.pop();
      vector<int> users = cover_domain[network_node];

      for (vector<int>::iterator it = users.begin(); it != users.end(); it++) {
        if (!pass[*it]) {
          pass[*it] = true;
          num--;

          if (num <= 0) {
            break;
          }

          for (vector<int>::iterator nit = server_candiate_locs[*it].begin();
               nit != server_candiate_locs[*it].end(); nit++) {
            connect_num[*nit]++;
            connect_num[*nit] += raoDong(  );

            Q.push(make_pair(connect_num[*nit], *nit));
          }
        }
      }
    }

    if (minServerNum > tempS1.locs.size()) {
      minServerNum = tempS1.locs.size();
    }

    bestLayoutFlow(tempS1);

    update(tempS1, true);

    server_candiate.push_back(tempS1);
  }
  
  
  if(value_lower> minServerNum * server_price ){
    value_lower = minServerNum * server_price;    
  }
  
}

void Loc_choose::supper_update() {

  
  
  for (size_t i = 0; i < server_candiate.size(); i++) {
    if (totCap == server_candiate[i].success_bw) {
      if (value_supper > server_candiate[i].total_price) {
        value_supper = server_candiate[i].total_price;
      }
    }
  }

  
}
void Loc_choose::forbidVirLinks() {
  for (int network_node = 0; network_node < network_node_num; network_node++) {
    graph.setWeight(2 * network_node, INF);
    graph.setWeight(2 * network_node + 1, INF);
  }
  for (int link = target_link_start; link < graph.getLink_num(); link++) {
    graph.setWeight(link, INF);
  }
}
void Loc_choose::releaseVirLinks() {
  for (int network_node = 0; network_node < network_node_num; network_node++) {
    graph.setWeight(2 * network_node, 0);
    graph.setWeight(2 * network_node + 1, 0);
  }
  for (int link = target_link_start; link < graph.getLink_num(); link++) {
    graph.setWeight(link, 0);
  }
}

bool Loc_choose::smallestUer() {
  if (2 + choosedServer.size() >= user_node_num) {
    Server_loc tempS;
    bestLayoutFlow(tempS);
    server_candiate.push_back(tempS);

    for (int user_node = 0; user_node < user_node_num; user_node++) {
      if (!user_direct_server[user_node]) {
        int network_node = user_to_network_map[user_node];
        Server_loc tempS;
        tempS.locs.insert(network_node);
        bestLayoutFlow(tempS);
        server_candiate.push_back(tempS);
      }
    }

    return true;
  }

  return false;
}
void Loc_choose::tryKServer(const int k) {}


void Loc_choose::bestLayoutFlow(Server_loc &server) {
  server.locs.insert(choosedServer.begin(), choosedServer.end());

  vector<int> caps = orignal_caps;
  for (set<int>::iterator nid = server.locs.begin(); nid != server.locs.end();
       nid++) {
    caps[*nid] = totCap;
  }

  vector<int> outs;
  vector<int> ins;
  vector<int> node_sum_value;
  server.user_in_bandwidth.resize(user_node_num, 0);
  pair<int, int> one_elem = graph.getMinCostMaxFlow(
      virtual_source, virtual_target, caps, outs, ins, node_sum_value);

  server.success_bw = one_elem.first;

  server.locs.clear();

  server.reach_node_value.clear();
  server.outBandwidth.clear();

  for (int network_node = 0; network_node < network_node_num; network_node++) {
    if (outs[network_node] > 0) {
      server.outBandwidth[network_node] = outs[network_node];
      server.locs.insert(network_node);
    }

    if (ins[network_node] > 0) {
      int user_node = network_node_user_map[network_node];
      server.user_in_bandwidth[user_node] = ins[network_node];
    }

    if (node_sum_value[network_node] > 0) {
      server.reach_node_value[network_node] = node_sum_value[network_node];
    }
  }
  server.total_price = one_elem.second + server.locs.size() * server_price;
}

void  Loc_choose::delete_canduate( void ){

  if(server_candiate.size(  )<2  ){
    return ;
  }
  
  vector<Server_loc> temps;
  Server_loc last=server_candiate[ 0 ];
  temps.push_back( last );
  for( size_t  i=1;i< server_candiate.size(  ); i++ ){
    if(last.locs!= server_candiate[ i ].locs  ){
      last=server_candiate[ i ];
      temps.push_back( last );
    }
  }
  server_candiate=temps;
  
}

void Loc_choose::update(Server_loc &server, bool recursive) {
  while (true) {
    int large_loc = -1;
    int large_value = -1;

    for (map<int, int>::iterator it = server.reach_node_value.begin();
         it != server.reach_node_value.end(); it++) {
      if (it->first > large_value) {
        large_value = it->first;
        large_loc = it->second;
      }
    }

    if (large_value >= server_price) {
      server.locs.insert(large_loc);
      bestLayoutFlow(server);
    } else {
      break;
    }
  }

  if (recursive) {
    bool state = true;
    while (state) {
      state = false;
      int small = INF;
      int smallNode = -1;
      for (map<int, int>::iterator it = server.outBandwidth.begin();
           it != server.outBandwidth.end(); it++) {
        if (it->second < small) {
          small = it->second;
          smallNode = it->first;
        }
      }

      if (small < smallest_user * delete_para) {
        Server_loc tempS = server;
        tempS.locs.erase(smallNode);
        bestLayoutFlow(tempS);
        update(tempS, false);

        if (tempS.success_bw > server.success_bw) {
          if (tempS.total_price <= server.total_price) {
            server = tempS;
            state = true;
          } else {
            server_candiate.push_back(tempS);
          }
        }

        if (tempS.success_bw == server.success_bw) {
          if (tempS.total_price < server.total_price) {
            server = tempS;
            state = true;
          }
        }
      }
    }
  }
}

/**
 *  @brief convert result to format of file
 *
 *
 * @return dynamic allow memory for output file
 */
char *Loc_choose::output() {
  int best_loc = -1;
  int best_value = INF;
  for (size_t i = 0; i < server_candiate.size(); i++) {
    if (totCap == server_candiate[i].success_bw) {
      if (best_value > server_candiate[i].total_price) {
        best_loc = i;
        best_value = server_candiate[i].total_price;
      }
    }
  }
  vector<int> caps = orignal_caps;
  for (set<int>::iterator it = server_candiate[best_loc].locs.begin();
       it != server_candiate[best_loc].locs.end(); it++) {
    caps[*it] = totCap;
  }

  vector<vector<int>> node_paths;
  vector<int> bws;

  graph.getMinCostMaxFlow(virtual_source, virtual_target, caps, node_paths,
                          bws);

  size_t len = 1024;
  for (size_t i = 0; i < node_paths.size(); i++) {
    len += 8 * node_paths[i].size();
  }
  len += 8 * node_paths.size();
  char *topo_file = new char[len];
  fill(topo_file, topo_file + len - 1, 0);
  sprintf(topo_file, "%d\n", node_paths.size());
  size_t start = strlen(topo_file);

  for (size_t i = 0; i < node_paths.size(); i++) {
    topo_file[start] = '\n';
    start++;
    size_t j = 1;
    int last = 0;
    for (size_t j = 1; j + 1 < node_paths[i].size(); j++) {
      last = node_paths[i][j];
      sprintf(topo_file + start, "%d ", last);
      start = strlen(topo_file);
    }

    sprintf(topo_file + start, "%d ", network_node_user_map[last]);
    start = strlen(topo_file);
    sprintf(topo_file + start, "%d", bws[i]);
    start = strlen(topo_file);
  }

  return topo_file;
}

void Loc_choose::initial_candidate_loc() {
  forbidVirLinks();

  for (int user_node = 0; user_node < user_node_num; user_node++) {
    server_candiate_locs[user_node].clear();

    int network_node = user_to_network_map[user_node];

    if (user_direct_server[user_node]) {
      server_candiate_locs[user_node].push_back(network_node);
      continue;
    }
    int limit=0;
    if(0!=user_demand[user_node]){
      limit = server_price / user_demand[user_node];
    }

    vector<pair<int, int>> tree;
    graph.dijkstra_limit_tree(network_node, limit, tree);

    for (vector<pair<int, int>>::iterator it = tree.begin(); it != tree.end();
         it++) {
      server_candiate_locs[user_node].push_back(it->first);
    }

    if (1 == server_candiate_locs[user_node].size()) {
      user_direct_server[user_node] = true;
      choosedServer.insert(network_node);
    }

    sort(server_candiate_locs[user_node].begin(),
         server_candiate_locs[user_node].end());
  }

  releaseVirLinks();
}

bool Loc_choose::domain_intersection_check() {
  
  allChoose.clear();
  bool re = true;


  for (int user_node = 0; user_node < user_node_num; user_node++) {
    std::vector<int> union_data;
    set_union(server_candiate_locs[user_node].begin(),
              server_candiate_locs[user_node].end(), allChoose.begin(),
              allChoose.end(), std::back_inserter(union_data));
    if (union_data.size() !=
        (server_candiate_locs[user_node].size() + allChoose.size())) {
      re = false;
    }
    sort(union_data.begin(), union_data.end());
    allChoose = union_data;
  }

  if (re) {
    return true;
  }

  vector<vector<int>> cover_domain(network_node_num);

  for (int user_node = 0; user_node < user_node_num; user_node++) {
    for (vector<int>::iterator it = server_candiate_locs[user_node].begin();
         it != server_candiate_locs[user_node].end(); it++) {
      cover_domain[*it].push_back(user_node);
    }
  }

  Fixed_heap Q( network_node_num);

  vector<int> connect_num(network_node_num);

  for (int network_node = 0; network_node < network_node_num; network_node++) {
    int len = cover_domain[network_node].size();
    Q.push(make_pair(-len, network_node));
    connect_num[network_node] = -len;
  }

  Server_loc tempS;

  vector<bool> pass(user_node_num, false);
  int num = user_node_num;

  while (!Q.empty()) {
    pair<int, int> p = Q.top();
    int network_node = p.second;

    Q.pop();
    vector<int> users = cover_domain[network_node];

    for (vector<int>::iterator it = users.begin(); it != users.end(); it++) {
      if (!pass[*it]) {

        pass[*it] = true;
        num--;

        tempS.locs.insert(network_node);
        if (num <= 0) {
          break;
        }

        for (vector<int>::iterator nit = server_candiate_locs[*it].begin();
             nit != server_candiate_locs[*it].end(); nit++) {
          connect_num[*nit]++;
          Q.push(make_pair(connect_num[*nit], *nit));
        }
      }
    }
  }

  int minServerNum = tempS.locs.size();

  bestLayoutFlow(tempS);

  update(tempS, true);

  server_candiate.push_back(tempS);

  for (int i = 0; i < (user_node_num) / 5 + 10; i++) {
    
    num = user_node_num;
    fill(pass.begin(  ), pass.end(  ), false);
    
    Q.clear(  );
    for (int network_node = 0; network_node < network_node_num;
         network_node++) {
      int len = cover_domain[network_node].size() +raoDong(  );

      Q.push(make_pair(-len, network_node));
      connect_num[network_node] = -len;
    }

    Server_loc tempS1;

    vector<bool> pass(user_node_num, false);
    int num = user_node_num;

    while (!Q.empty()) {
      pair<int,int> p = Q.top();
      int network_node = p.second;


      Q.pop();
      vector<int> users = cover_domain[network_node];

      for (vector<int>::iterator it = users.begin(); it != users.end(); it++) {
        if (!pass[*it]) {

          pass[*it] = true;
          num--;

          tempS1.locs.insert(network_node);
          if (num <= 0) {
            break;
          }

          for (vector<int>::iterator nit = server_candiate_locs[*it].begin();
               nit != server_candiate_locs[*it].end(); nit++) {
            connect_num[*nit]++;
            connect_num[*nit] += raoDong(  );

            Q.push(make_pair(connect_num[*nit], *nit));
          }
        }
      }
    }

    if (minServerNum > tempS1.locs.size()) {
      minServerNum = tempS1.locs.size();
    }

    bestLayoutFlow(tempS1);

    update(tempS1, true);

    server_candiate.push_back(tempS1);
  }
  
  if(value_lower> minServerNum * server_price ){
    value_lower = minServerNum * server_price;    
  }


  return false;
}
char *Loc_choose::solve() {
  /**
   *  there is no user node
   *
   */
  value_lower=INF;
  value_supper=0;
  
  if (0 == user_node_num) {
    char *topo_file = new char[1024];
    fill(topo_file, topo_file + 1023, 0);
    sprintf(topo_file, "NA");
    return topo_file;
  }

  server_candiate.clear();

  initial_candidate_loc();

  Server_loc tempS;

  for (int user_node = 0; user_node < user_node_num; user_node++) {
    int network_node = user_to_network_map[user_node];
    tempS.locs.insert(network_node);
  }

  bestLayoutFlow(tempS);

  server_candiate.push_back(tempS);

  if (domain_intersection_check()) {
    return output();
  }

  if (smallestUer()) {
    return output();
  }

  value_supper = user_node_num * server_price;

  for (vector<int>::iterator it = allChoose.begin(); it != allChoose.end();
       it++) {
    int network_node = *it;

    if (choosedServer.find(network_node) != choosedServer.end()) {
      continue;
    }

    Server_loc temp;
    temp.locs.insert(network_node);
    bestLayoutFlow(temp);
    if ((totCap == temp.success_bw) && (temp.total_price < value_supper)) {
      value_supper = temp.total_price;
    }
    server_candiate.push_back(temp);
  }

  /**
   *  one server can support
   *
   */
  if (value_supper <= 2 * server_price) {
    return output();
  }

  if (value_lower < 2 * server_price) {
    value_lower = 2 * server_price;
  }

  for (int k = 0; k < user_node_num/10+5; k++) {

    sort(server_candiate.begin(), server_candiate.end());

    delete_canduate(  );
    
    lower_update();
    
    supper_update();
    
    if ((value_supper / server_price) == (value_lower / server_price)) {
      return output();
    }
    
  }

  return output();
}
}
