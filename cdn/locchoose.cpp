
#include<iostream>

#include<string.h>

#include<algorithm>

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
  double re=clock();
  return re/CLOCKS_PER_SEC;
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
      }else{
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
  for(int network_node=0; network_node< network_node_num; network_node++){
      vector<int> dis; 
      graph.dijkstra_tree(network_node, dis);
      for(int user_node=0; user_node< user_node_num; user_node++){
          int temp=user_to_network_map[user_node];
          network_to_user_inv_distance[network_node][user_node]=1.0/(dis[temp]+1); 
      }
  }
  releaseVirLinks();
}

void Loc_choose::initial_case(void){
  vector<vector<int>> network_node_cover_domain(network_node_num);

  for (int user_node = 0; user_node < user_node_num; user_node++) {
    for (vector<int>::iterator it = server_candiate_locs[user_node].begin();
         it != server_candiate_locs[user_node].end(); it++) {
      network_node_cover_domain[*it].push_back(user_node);
    }
  }

  Fixed_heap Q( network_node_num);

  vector<int> connect_num(network_node_num);

  for (int network_node = 0; network_node < network_node_num; network_node++) {
    int len = network_node_cover_domain[network_node].size();
    Q.push(make_pair(-len, network_node));
    connect_num[network_node] = -len;
  }

  
  vector<bool> pass(user_node_num, false);
  int num = user_node_num;
  
  int minServerNum =  INF;
  
  for (int i = 0; i< 20; i++) {
    
    num = user_node_num;
    fill(pass.begin(  ), pass.end(  ), false);
    
    Q.clear(  );
    for (int network_node = 0; network_node < network_node_num;
         network_node++) {
      int len = network_node_cover_domain[network_node].size() +raoDong(  );

      Q.push(make_pair(-len, network_node));
      connect_num[network_node] = -len;
    }

    Server_loc tempServer;

    int num = user_node_num;

    while (!Q.empty()) {
      pair<int,int> p = Q.top();
      int network_node = p.second;
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
            connect_num[*nit] += raoDong(  );

            Q.push(make_pair(connect_num[*nit], *nit));
          }
        }
      }
    }

    if (minServerNum > tempServer.locs.size()) {
      minServerNum = tempServer.locs.size();
    }

    bestLayoutFlow(tempServer);

    update(tempServer, true);

    server_candiate.push_back(tempServer);
  }
  
  
  if(value_lower> minServerNum * server_price ){
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
      temp_value += server_candiate[i].total_price *
                    ((totCap - cap) / (server_candiate[i].success_bw + 0.01));
      if (value_lower > temp_value) {
        value_lower = temp_value;
      }
      break;
    }
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
      if(temps.size()>= 40){
        break;
      }
    }
  }
  server_candiate=temps;
}

 void Loc_choose::generate_case( Server_loc& lhs, Server_loc &rhs){
  
    set<int> locs=lhs.locs;
    
    int len1=locs.size();
        
    int len2=rhs.locs.size();
      
    locs.insert(rhs.locs.begin(), rhs.locs.end());
  
    vector<int> tempLoc(locs.begin(), locs.end());
    random_shuffle(tempLoc.begin(), tempLoc.end());
    Server_loc tempS;
    int minL= len1< len2? len1:len2;
    minL--; ;
    if(minL<0) {
      minL=0;
    }

    int maxL= len1< len2? len2:len1;
    maxL++;
    if(maxL>tempLoc.size()){
      maxL=tempLoc.size();
    }

    int outLen = minL + (rand() % (int)(maxL - minL + 1));
      
    tempS.locs.insert(tempLoc.begin(), tempLoc.begin()+outLen);
    bestLayoutFlow(tempS);
    if(network_node_num*user_node_num<100000){
        
        update(tempS, deleteSmall);
    }
    server_candiate.push_back(tempS);
  
}
void  Loc_choose::randdom_generate(void){

  size_t firstLen=first_class_candiate_num;
  
  if(server_candiate.size()< firstLen){
    firstLen=server_candiate.size();
  }
  
  for(size_t i=0; i< firstLen; i++){

    for(size_t j=i; j< firstLen; j++){

      generate_case(server_candiate[i], server_candiate[j]);
    }
  }

  sort(server_candiate.begin(), server_candiate.end());
    
  int secondLen=server_candiate.size();
  if(secondLen>0){
      secondLen--;
  }

  
  for(size_t i=0; i< firstLen; i++){

    int rid = firstLen + (rand() % (int)(secondLen - firstLen + 1));
    generate_case(server_candiate[i], server_candiate[rid]);

    rid = firstLen + (rand() % (int)(secondLen - firstLen + 1));
      
    generate_case(server_candiate[i], server_candiate[rid]);

  }
}

  void Loc_choose::update(Server_loc &server, bool recursive) {

    while(server.success_bw< totCap){
      double left=totCap-server.success_bw;
      int es_add_num=1;//((left/(totCap-server.success_bw+0.1))*server.locs.size())/5+1;
  
      vector<pair<double, int> > candiateN;
        
      vector<int> leftBandwidth(user_node_num, 0);
      for(int user_node=0; user_node< user_node_num;  user_node++){
        leftBandwidth[user_node]=user_demand[user_node]-server.user_in_bandwidth[user_node];
      }

      for(int network_node=0; network_node< network_node_num; network_node++){
        if(server.locs.find(network_node)==server.locs.end()){
          float temp_value=0;
          for(int user_node=0; user_node< user_node_num; user_node++){
            temp_value+=leftBandwidth[user_node]*network_to_user_inv_distance[network_node][user_node];
          }
          candiateN.push_back(make_pair(temp_value, network_node));
        }
      }
    
      sort(candiateN.rbegin(), candiateN.rend());
      for(size_t i=0; i< candiateN.size() && i< es_add_num; i++){
        server.locs.insert(candiateN[i].second);      
      }
    
      bestLayoutFlow(server);
    
    }
  
  
    while (true) {

      vector<pair<int, int> > candiateN;
      for (map<int, int>::iterator it = server.reach_node_value.begin();
           it != server.reach_node_value.end(); it++) {
        if ( it->second >=server_price) {
          candiateN.push_back(make_pair(it->second, it->first));
        }
      }

      if(!candiateN.empty()){
        sort(candiateN.rbegin(), candiateN.rend());
        int addN=1;//candiateN.size()/5+1;
        for(size_t i=0; i< addN && i< candiateN.size(); i++){
          server.locs.insert(candiateN[i].second);
        }

        bestLayoutFlow(server);

      }else{
        break;
      }
    }

    if (recursive) {
      bool state=true;

      while (state) {
      
        state=false;
        vector<pair<int, int> > candiateN;
        for (map<int, int>::iterator it = server.outBandwidth.begin();
             it != server.outBandwidth.end(); it++) {
          if (it->second  <= smallest_user * delete_para) {
            candiateN.push_back(make_pair(it->second, it->first ));
          }
        }
      
        if(!candiateN.empty()){
        
          sort(candiateN.begin(), candiateN.end());
        
          int deleteN=1;//candiateN.size()/5+1;
        
          Server_loc tempS = server;
        
          for(size_t i=0; i< deleteN && i< candiateN.size(); i++){
            tempS.locs.erase(candiateN[i].second);
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
  if(best_loc>-1 && server_candiate[best_loc].locs.empty()){
    char *topo_file = new char[1024];
    fill(topo_file, topo_file + 1023, 0);
    sprintf(topo_file, "NA");
    return topo_file;
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

  vector<int> choose_time(network_node_num, 0);
  
  for (int user_node = 0; user_node < user_node_num; user_node++) {

    for(vector<int>::iterator it=server_candiate_locs[user_node].begin(); it!=server_candiate_locs[user_node].end(); it++){
      choose_time[*it]++;
    }
  }
  
  bool re = true;
  for(int network_node=0; network_node< network_node_num; network_node++){
    if(choose_time[network_node]>0){
      allChoose.push_back(network_node);
      if(choose_time[network_node]>1){
        re=false;
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

  Fixed_heap Q( network_node_num);

  vector<int> connect_num(network_node_num);

  for (int network_node = 0; network_node < network_node_num; network_node++) {
    int len = network_node_cover_domain[network_node].size();
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
          Q.push(make_pair(connect_num[*nit], *nit));
        }
      }
    }
  }

  int minServerNum = tempS.locs.size();

  bestLayoutFlow(tempS);

  update(tempS, true);

  server_candiate.push_back(tempS);
  
  
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
  
  if ((0 == user_node_num)|| (0==totCap) ) {
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
  server_candiate.push_back(tempS);
  if(0==server_price){
    return output();
  }

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

  int steable_time=0;
  int last_best_value=-1;
  for (int k = 0; k < user_node_num+5; k++) {
    std::cout << "time(s): "<<systemTime()-start_time<<" value: "<<value_supper << std::endl;
    if(systemTime()-start_time> time_bound-10){
      break;
    }
    // if(systemTime()-start_time> 300){
    //   break;
    // }
    if(value_supper==last_best_value){
      steable_time++;
      first_class_candiate_num--;
      if(first_class_candiate_num<1){
         first_class_candiate_num=1;
      }
    }else{
      deleteSmall=false;
      steable_time=0;
    }
    if(steable_time>5){
        first_class_candiate_num=3*steable_time/5; 
    }
    
    if(steable_time>10){

       deleteSmall=true;
    }

            
    if(steable_time>20){
      break;
    }
    
    last_best_value=value_supper;

    
    sort(server_candiate.begin(), server_candiate.end());

    delete_canduate(  );

    randdom_generate();
    
    supper_update();
    
    if ((value_supper / server_price) == (value_lower / server_price)) {
      return output();
    }
    
  }

  return output();
}
}
