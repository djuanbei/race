
#include <iostream>

#include <string.h>

#include <algorithm>

#include <assert.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include "locchoose.h"

namespace raptor {
using namespace std;

double systemTime(void) {
  struct timespec start;
#ifndef __MACH__
  clock_gettime(CLOCK_MONOTONIC, &start);
#else
  double re = clock();
  return re / CLOCKS_PER_SEC;
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

      if (sumBw < 2 * user_demand[user_node]) {
        user_direct_server[user_node] = true;
        choosedServer.insert(network_node);
        tempC -= user_demand[user_node];
      } else {
        demandC.push_back(make_pair(user_demand[user_node], user_node));
      }
    }
  }

  sort(demandC.rbegin(), demandC.rend());
  for (size_t i = 0; i + 1 < demandC.size(); i++) {
    /**
     *  when one user demand greater or equal then sum of all
     * other users demand then this user direct connect node must be
     * a server
     *
     */
    int user_node = demandC[i].second;  // it->second;
    int network_node = user_to_network_map[user_node];

    if (!user_direct_server[user_node]) {
      if (2 * (demandC[i].first) >= tempC) {
        user_direct_server[user_node] = true;
        choosedServer.insert(network_node);
        tempC -= demandC[i].first;
      } else {
        break;
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

  forbidVirLinks();
  for (int network_node = 0; network_node < network_node_num; network_node++) {
    vector<int> dis;
    graph.dijkstra_tree(network_node, dis);
    for (int user_node = 0; user_node < user_node_num; user_node++) {
      int temp = user_to_network_map[user_node];
      network_to_user_inv_distance[network_node][user_node] =
        1.0 / (dis[temp] + 1);
    }
  }
  releaseVirLinks();
}

void Loc_choose::initial_case(void) {
  vector<vector<int>> network_node_cover_domain(network_node_num);

  for (int user_node = 0; user_node < user_node_num; user_node++) {
    for (vector<int>::iterator it = server_candiate_locs[user_node].begin();
         it != server_candiate_locs[user_node].end(); it++) {
      network_node_cover_domain[*it].push_back(user_node);
    }
  }

  Fixed_heap Q(network_node_num);

  vector<int> connect_num(network_node_num);

  for (int network_node = 0; network_node < network_node_num; network_node++) {
    int len = network_node_cover_domain[network_node].size();
    Q.push(-len, network_node);
    connect_num[network_node] = -len;
  }

  vector<bool> pass(user_node_num, false);
  int num = user_node_num;

  int minServerNum = INF;

  for (int i = 0; i < para.initcase_num; i++) {
    num = user_node_num;
    fill(pass.begin(), pass.end(), false);

    Q.clear();
    for (int network_node = 0; network_node < network_node_num;
         network_node++) {
      int len = network_node_cover_domain[network_node].size() + raoDong();

      Q.push(-len, network_node);
      connect_num[network_node] = -len;
    }

    Server_loc tempServer;

    int num = user_node_num;

    while (!Q.empty()) {
      int temp;
      int network_node;
      // pair<int, int> p =
      Q.top(temp, network_node);
      // int network_node = p.second;
      tempServer.locs.insert(network_node);

      Q.pop();
      vector<int> users = network_node_cover_domain[network_node];

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
            connect_num[*nit] += raoDong();

            Q.push(connect_num[*nit], *nit);
          }
        }
      }
    }

    if (minServerNum > tempServer.locs.size()) {
      minServerNum = tempServer.locs.size();
    }

    bestLayoutFlow(tempServer);

    update(tempServer, true);
    addLoc(tempServer);
  }

  if (value_lower > minServerNum * server_price) {
    value_lower = minServerNum * server_price;
  }
}
void Loc_choose::lower_update() {
  sort(server_candiate.begin(), server_candiate.end());

  int cap = 0;
  int temp_value = 0;
  for (size_t i = 0; i < server_candiate.size(); i++) {
    if (cap + server_candiate[i].success_bw < totCap) {
      cap += server_candiate[i].success_bw;
      temp_value += server_candiate[i].total_price;
    } else {
      temp_value +=
          server_candiate[i].total_price *
          ((totCap - cap) / (server_candiate[i].success_bw + 0.01));
      if (value_lower > temp_value) {
        value_lower = temp_value;
      }
      break;
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
    addLoc(tempS);

    for (int user_node = 0; user_node < user_node_num; user_node++) {
      if (!user_direct_server[user_node]) {
        int network_node = user_to_network_map[user_node];
        Server_loc tempS;
        tempS.locs.insert(network_node);
        bestLayoutFlow(tempS);
        addLoc(tempS);
      }
    }

    return true;
  }

  return false;
}

 bool Loc_choose::bestLayoutFlow(Server_loc &server, bool more_check) {
  server.locs.insert(choosedServer.begin(), choosedServer.end());

  for( size_t i=0; i< server_candiate.size(  ); i++ ){
    if( server.locs==server_candiate[ i ].locs){
      server=server_candiate[ i ];
      return time_end(  );
    }
  }

  vector<int> caps = orignal_caps;
  int leftCap = totCap;
  for (set<int>::iterator nid = server.locs.begin(); nid != server.locs.end();
       nid++) {
    caps[*nid] = totCap;
    fill( sum_of_pass_flow[ *nid ].begin(  ) ,sum_of_pass_flow[ *nid].end(  ), 0);
    
    if (network_node_user_map[*nid] > -1) {
      int user_node = network_node_user_map[*nid];
      int inlink = target_link_start / 2 + user_node;
      leftCap -= user_demand[user_node];

      caps[inlink] = 0;
    }
  }

  vector<int> outs;
  vector<int> ins;
  vector<int> node_sum_value;
  server.user_in_bandwidth.resize(user_node_num, 0);

  pair<int, int> one_elem = graph.getMinCostMaxFlow(
      virtual_source, virtual_target, caps, leftCap, outs, ins, node_sum_value, sum_of_pass_flow );

  server.success_bw = one_elem.first;

  server.reach_node_value.clear();
  server.outBandwidth.clear();

  for (set<int>::iterator nid = server.locs.begin(); nid != server.locs.end();
       nid++) {
    if (network_node_user_map[*nid] > -1) {
      int user_node = network_node_user_map[*nid];
      ins[*nid] += user_demand[user_node];
      assert(ins[*nid] == user_demand[user_node]);

      server.success_bw += user_demand[user_node];

      outs[*nid] += user_demand[user_node];
    }
  }

  server.locs.clear();

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

    
  if( !time_end()&& totCap==server.success_bw ){

    set<int> locs=server.locs;
    for(set<int>::iterator it= locs.begin(  ); it!= locs.end(  ); it++ ){
      Server_loc tempS=server;
      int sum=server.outBandwidth[*it];
      int node=-1;
      int big=-1;
      
      for(int network_node =0; network_node< network_node_num; network_node++){
        if(2*sum_of_pass_flow[ *it ][ network_node ]>=sum  ){
          if(sum_of_pass_flow[ *it ][ network_node ]>big){
            big=sum_of_pass_flow[ *it ][ network_node ];
            node=network_node;
          }
        }
      }
      
      if(node>-1){
        server.locs.erase( *it );
        server.locs.insert( node );
        addLoc(tempS);

        return bestLayoutFlow( server );
      }
    }
  }
  if(more_check){
    if( !time_end()&& totCap==server.success_bw ){

      Server_loc tempS=server;

      set<int> locs=server.locs;
      for(set<int>::iterator it= locs.begin(  ); it!= locs.end(  ); it++ ){
        int sum=server.outBandwidth[*it];
        int node=-1;
        int big=-1;
      
        for(int network_node =0; network_node< network_node_num; network_node++){
          if(sum_of_pass_flow[ *it ][ network_node ]>big){
            big=sum_of_pass_flow[ *it ][ network_node ];
            node=network_node;
          }

        }
      
        if(node>-1 &&3*big>sum ){
          Server_loc temp=server;
          temp.locs.erase( *it );
          temp.locs.insert( node );
          addLoc(tempS);
          bestLayoutFlow( temp, false );
          if(time_end()){
            return true;
          }
                
        }
      }

    }
  }
  

  return time_end();
}
  


void Loc_choose::delete_candiate(void) {

  deleteRepeat(server_candiate, para.part_left_num);

}

void Loc_choose::generate_case(Server_loc &lhs, Server_loc &rhs) {
  set<int> locs = lhs.locs;

  int len1 = locs.size();

  int len2 = rhs.locs.size();

  locs.insert(rhs.locs.begin(), rhs.locs.end());

  bool same = (locs.size() == len1);

  vector<int> tempLoc(locs.begin(), locs.end());
  random_shuffle(tempLoc.begin(), tempLoc.end());

  int minL = len1 < len2 ? len1 : len2;
  minL--;

  if (minL < 0) {
    minL = 0;
  }

  int maxL = len1 < len2 ? len2 : len1;
  maxL ++;
  if (maxL > tempLoc.size()) {
    maxL = tempLoc.size();
  }

  int outLen = minL + (rand() % (int)(maxL - minL + 1));
  if (same) {
    if (outLen == tempLoc.size()) {
      outLen--;
    }
    if (outLen <= 0) {
      return;
    }
  }
  set<int> locSet(tempLoc.begin(), tempLoc.begin() + outLen  );
  
  bool add=true;
  for( size_t i= 0; i<randCase.size(  ); i++ ){
    if( locSet==randCase[ i ] ){
      add=false;
      break;
    }
  }
  
  if( add ){
    randCase.push_back( locSet );
  }
}



void Loc_choose::randdom_generate(void) {
  randCase.clear(  );
  size_t firstLen = para.first_class_candiate_num;
  sort(server_candiate.begin(), server_candiate.end());

  if (server_candiate.size() < firstLen) {
    firstLen = server_candiate.size();
  }

  for (size_t i = 0; i < firstLen; i++) {
    if(para.randAddNum>0){
      Server_loc tempS=server_candiate[ i ];
      for( int k=0; k<para.randAddNum ; k++  ){

        int network_node =rand(  )% network_node_num;
        while(server_candiate[ i ].locs.find(network_node  )!= server_candiate[ i ].locs.end(  ) )      {
          network_node =rand(  )% network_node_num;
        }
        tempS.locs.insert( network_node );
      }
      
      bestLayoutFlow( tempS );
      update( tempS, true );
      addLoc( tempS );
    }

      
    for (size_t j = i; j < firstLen; j++) {
      generate_case(server_candiate[i], server_candiate[j]);
    }
  
  }
  
  for(size_t i=0; i+1< best_loc.size(); i++){
    int best=-1;
    double value=2;
    for(size_t j=i+1; j< best_loc.size(); j++){
      int temp=simliar1(best_loc[i], best_loc[j]);
      if(temp<value){
        value=temp;
        best=i;
      }
    }
    
    if(best>-1){
      generate_case(best_loc[i], best_loc[best]);
    }
  }
  
  sort(server_candiate.begin(), server_candiate.end());
  
  if (server_candiate.size() > 0&& best_loc.size()>0) {
    int secondLen = server_candiate.size()/2;

    secondLen--;
    
    int minL=firstLen;
    firstLen=0;
    for (size_t i = 0; i < minL; i++) {

      int rid = firstLen + (rand() % (int)(secondLen - firstLen + 1));
      generate_case(best_loc[0], server_candiate[rid]);

      rid = firstLen + (rand() % (int)(secondLen - firstLen + 1));
      generate_case(best_loc[0], server_candiate[rid]);

    }
  }
  
  
  for( size_t i=0; i< randCase.size(  ); i++ ){
    Server_loc tempS;
    tempS.locs=randCase[ i ];
    bestLayoutFlow(tempS);
    if (network_node_num * user_node_num < 100000) {
      update(tempS, false);
    } else {
      update(tempS, para.deleteSmall);
    }
    addLoc(tempS);
    if(time_end()){
       return;
    }
  }
}

void Loc_choose::update(Server_loc &server, bool recursive) {

  while (server.success_bw < totCap) {
    
    // double left = totCap - server.success_bw;
    
    int es_add_num = 1;
    if (user_node_num * network_node_num >= para.large_scale) {
      es_add_num = para.add_num;
    }

    vector<pair<double, int>> candiateN;

    vector<int> leftBandwidth(user_node_num, 0);
    for (int user_node = 0; user_node < user_node_num; user_node++) {
      leftBandwidth[user_node] =
          user_demand[user_node] - server.user_in_bandwidth[user_node];
    }

    for (int network_node = 0; network_node < network_node_num;
         network_node++) {
      if (server.locs.find(network_node) == server.locs.end()) {
        float temp_value = 0;
        for (int user_node = 0; user_node < user_node_num; user_node++) {
          temp_value +=
              (leftBandwidth[user_node] *
               network_to_user_inv_distance[network_node][user_node]) *
              (1 + ((rand() % 100) / 1000.0));
        }
        candiateN.push_back(make_pair(temp_value, network_node));
      }
    }

    sort(candiateN.rbegin(), candiateN.rend());
    for (size_t i = 0; i < candiateN.size() && i < es_add_num; i++) {
      server.locs.insert(candiateN[i].second);
    }

    if (bestLayoutFlow(server)) {
      return;
    }
  }


  while (true) {
    vector<pair<int, int>> candiateN;
    for (map<int, int>::iterator it = server.reach_node_value.begin();
         it != server.reach_node_value.end(); it++) {
      if (it->second >= server_price) {
        candiateN.push_back(make_pair(it->second, it->first));
      }
    }
    if (!candiateN.empty()) {
      if(0==para.iterator_num){
        for(size_t i=0; i<candiateN.size(); i++){
          server.locs.insert(candiateN[i].second);
        }
        if (bestLayoutFlow(server)) {
          return;
        }
      }else{

        sort(candiateN.rbegin(), candiateN.rend());
        server.locs.insert(candiateN[0].second);
      
        if (bestLayoutFlow(server)) {
          return;
        }
      } 
    }else {
      break;
    }
  }

  if (recursive) {
    bool state = true;

    while (state) {
      if(time_end()){
        return;
      }
      
      state = false;
      vector<pair<int, int>> candiateN;
      for (map<int, int>::iterator it = server.outBandwidth.begin();
           it != server.outBandwidth.end(); it++) {
        if (it->second <= smallest_user * para.delete_para) {
          candiateN.push_back(make_pair(it->second, it->first));
        }
      }

      if (!candiateN.empty()) {
        sort(candiateN.begin(), candiateN.end());

        int deleteN = 1;

        Server_loc tempS = server;

        for (size_t i = 0; i < deleteN && i < candiateN.size(); i++) {
          tempS.locs.erase(candiateN[i].second);

          if (bestLayoutFlow(tempS)) {
            return;
          }
          update(tempS, false);

          if (tempS.success_bw > server.success_bw) {
            if (tempS.total_price <= server.total_price) {
              server = tempS;
              state = true;
            } else {
              addLoc(tempS);
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
}

/**
 *  @brief convert result to format of file
 *
 *
 * @return dynamic allow memory for output file
 */
char *Loc_choose::output() {
  
  if (best_loc.empty()) {
    char *topo_file = new char[1024];
    fill(topo_file, topo_file + 1023, 0);
    sprintf(topo_file, "NA");
    return topo_file;
  }
  vector<int> caps = orignal_caps;
  for (set<int>::iterator it = best_loc.front().locs.begin(); it != best_loc.front().locs.end();
       it++) {
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
    int limit = 0;
    if (0 != user_demand[user_node]) {
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

  vector<int> choose_time(network_node_num, 0);

  for (int user_node = 0; user_node < user_node_num; user_node++) {
    for (vector<int>::iterator it = server_candiate_locs[user_node].begin();
         it != server_candiate_locs[user_node].end(); it++) {
      choose_time[*it]++;
    }
  }

  bool re = true;
  for (int network_node = 0; network_node < network_node_num; network_node++) {
    if (choose_time[network_node] > 0) {
      allChoose.push_back(network_node);
      if (choose_time[network_node] > 1) {
        re = false;
      }
    }
  }

  if (re) {
    return true;
  }

  vector<vector<int>> network_node_cover_domain(network_node_num);

  for (int user_node = 0; user_node < user_node_num; user_node++) {
    for (vector<int>::iterator it = server_candiate_locs[user_node].begin();
         it != server_candiate_locs[user_node].end(); it++) {
      network_node_cover_domain[*it].push_back(user_node);
    }
  }

  Fixed_heap Q(network_node_num);

  vector<int> connect_num(network_node_num);

  for (int network_node = 0; network_node < network_node_num; network_node++) {
    int len = network_node_cover_domain[network_node].size();
    Q.push(-len, network_node);
    connect_num[network_node] = -len;
  }

  Server_loc tempS;

  vector<bool> pass(user_node_num, false);
  int num = user_node_num;

  while (!Q.empty()) {
    int temp, network_node;
    // pair<int, int> p =
    Q.top(temp, network_node);
    // int network_node = p.second;

    Q.pop();
    vector<int> users = network_node_cover_domain[network_node];

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
          Q.push(connect_num[*nit], *nit);
        }
      }
    }
  }

  int minServerNum = tempS.locs.size();

  bestLayoutFlow(tempS);

  update(tempS, true);
  addLoc(tempS);


  for (vector<int>::iterator it = allChoose.begin(); it != allChoose.end();
       it++) {
    int network_node = *it;

    if (choosedServer.find(network_node) != choosedServer.end()) {
      continue;
    }
    
    if(network_node_user_map[*it]>-1){
      continue;
    }

    Server_loc temp;
    temp.locs.insert(network_node);
    bestLayoutFlow(temp);

    if ((totCap == temp.success_bw) && (temp.total_price < value_supper)) {
      value_supper = temp.total_price;
    }
    addLoc(temp);
  }


  if (value_lower > minServerNum * server_price) {
    value_lower = minServerNum * server_price;
  }

  return false;
}
char *Loc_choose::solve() {
  /**
   *  there is no user node
   *
   */
  value_lower = INF;
  value_supper = 0;

  if ((0 == user_node_num) || (0 == totCap)) {
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
  update(tempS, true);
  addLoc(tempS);

  if (0 == server_price) {
    return output();
  }

  if (domain_intersection_check()) {
    return output();
  }

  if (smallestUer()) {
    return output();
  }

  value_supper = user_node_num * server_price;

  lower_update();

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
  initial_case();

  int steable_time = 0;
  int last_best_value = -1;
  while (true) {
    para.iterator_num++;
    
    std::cout << "time(s): " << systemTime() - start_time
              << " value: " << value_supper<<" server_value: "<<best_loc.front().locs.size()*server_price << std::endl;

    if (time_end()) {
      break;
    }
    if (network_node_num < 300) {
      initial_case();
    }

    if (value_supper == last_best_value) {
      para.randAddNum++;

      if( network_node_num<600 ){
        if( para.randAddNum>2 ){
          para.randAddNum=2;          
        }
      }
      
      if(para.randAddNum>value_supper/server_price){
        para.randAddNum=value_supper/server_price;
      }
      para.add_num--;
      if (para.add_num < 1) {
        para.add_num = 1;
      }
      steable_time++;
      para.first_class_candiate_num--;
      if (para.first_class_candiate_num < 1) {
        para.first_class_candiate_num = 1;
      }
    } else {
      para.randAddNum=0;
      para.add_num++;
      if (para.add_num > 5) {
        para.add_num = 5;
      }
      para.deleteSmall = false;
      steable_time = 0;
    }


    if (steable_time > 5) {

      para.first_class_candiate_num = 2 * steable_time / 5 + 1;
    }

    if (steable_time > 10) {

      para.deleteSmall = true;
    }

    if (steable_time > para.stable_bound) {
      break;
    }

    last_best_value = value_supper;

    delete_candiate();

    randdom_generate();
    

    if ((value_supper / server_price) == (value_lower / server_price)) {
      return output();
    }
  }

  return output();
}
}
