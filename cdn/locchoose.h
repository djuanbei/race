#ifndef _LOC_CHOOSE_H
#define _LOC_CHOOSE_H

#include <cstdlib>
#include <limits>
#include <map>
#include <set>
#include <utility>
#include<cmath>
#include <vector>
#include "graph.h"

namespace raptor {

using namespace std;

double systemTime(void);

const static int time_bound = 90;

class Loc_choose {
  struct Server_loc {
    int success_bw;
    int total_price;

    set<int> locs;
    map<int, int> outBandwidth;
    vector<int> user_in_bandwidth;

    map<int, int> reach_node_value;
    bool isCheck;
    Server_loc() : success_bw(0), total_price(0),isCheck(false)  {}
    bool operator<(const Server_loc &other) const {
      if (0 == success_bw) {
        return 0 == other.success_bw;
      }
      if (0 == other.success_bw) {
        return true;
      }
      return (total_price /(success_bw+0.1)) <
        (other.total_price /(other.success_bw +0.1));
    }
  };

  float simliar(const Server_loc &lhs, const Server_loc &rhs) const {
    double sum = 0;
    double fen1, fen2;
    fen1 = fen2 = 0;

    for (int i = 0; i < user_node_num; i++) {
      sum += lhs.user_in_bandwidth[i] * rhs.user_in_bandwidth[i];
      fen1 += lhs.user_in_bandwidth[i] * lhs.user_in_bandwidth[i];
      fen2 += rhs.user_in_bandwidth[i] * rhs.user_in_bandwidth[i];
    }
    fen1 = sqrt(fen1);
    fen2 = sqrt(fen2);
    return (sum) / (fen1 * fen2 + 0.1);
  }

  float simliar1(Server_loc &lhs, Server_loc &rhs) {
    double sum = 0;
    double fen1, fen2;
    fen1 = fen2 = 0;

    for (int i = 0; i < user_node_num; i++) {
      int network_node = user_to_network_map[i];

      sum += lhs.reach_node_value[network_node] *
             rhs.reach_node_value[network_node];
      fen1 += lhs.reach_node_value[network_node] *
              lhs.reach_node_value[network_node];
      fen2 += rhs.reach_node_value[network_node] *
              rhs.reach_node_value[network_node];
    }
    fen1 = sqrt(fen1);
    fen2 = sqrt(fen2);
    return (sum) / (fen1 * fen2 + 0.1); 
  }

  struct Para {
    int randTryNum;

    double delete_para;

    int first_class_candiate_num;

    bool deleteSmall;

    int stable_bound;
    int success_left_num;
    int part_left_num;
    int add_num;
    int large_scale;
    int initcase_num;
    int randAddNum;
    int iterator_num;
    Para() {
      randTryNum = 100;
      delete_para = 0.9;
      first_class_candiate_num = 1;
      deleteSmall = false;
      stable_bound = 10000;
      success_left_num = 3;
      part_left_num = 40;
      add_num = 7;
      large_scale = 150000;
      initcase_num = 30;
      randAddNum=0;
      iterator_num=0;
    }
  };

 private:
  undirected_graph &graph;

  int network_node_num;

  int user_node_num;

  int server_price;

  int virtual_source, virtual_target;

  int target_link_start;

  const vector<int> &network_node_user_map;
  const vector<int> &user_demand;

  int totCap;
  const vector<int> &orignal_caps;

  int value_supper, value_lower;

  vector<Server_loc> success_server_candiate;

  vector<Server_loc> server_candiate;

  int lest_link_price;

  int large_link_price;
  double mean_link_price;
  int middle_link_price;

  int smallest_user;
  int large_user;
  
  vector<int> choose_time;

  /**
   * for every user u there at less a network node in server_candiate_locs[u]
   * used as a server
   *
   */

  vector<vector<int>> server_candiate_locs;

  vector<int> user_to_network_map;

  vector<bool> user_direct_server;

  set<int> choosedServer;

  vector<int> allChoose;
  
  vector<set<int>> randCase;

  vector<vector<float>> network_to_user_inv_distance;

  vector<vector<int>> sum_of_pass_flow;
  
  Para para;

  void deleteRepeat(vector<Server_loc> &loc, const int num ){
    
    if(loc.size()<2){
      return;
    }
    
    sort(loc.begin(), loc.end());
    set<int> last=loc.front().locs;
    vector<Server_loc> temp;
    temp.push_back(loc.front());
    
    for(size_t i=1; i<loc.size(); i++ ){
      if(last!=loc[i].locs){
        last=loc[i].locs;
        temp.push_back(loc[i]);
        if(temp.size()>=num){
          break;
        }
      }
    }
    
    loc=temp;
  }
  
  void addLoc(Server_loc &loc) {
    
    if (totCap == loc.success_bw) {
      success_server_candiate.push_back(loc);
      if(loc.total_price<best_loc.total_price){
        value_supper=best_loc.total_price;
        best_loc=loc;
      }
    }
    server_candiate.push_back(loc);
    for( set<int>::iterator it=loc.locs.begin(  ); it!= loc.locs.end(  ); it++ ){
      choose_time[ *it ]++;
    }

  }

  bool time_end() const { return systemTime() - start_time > time_bound -3; }

  int raoDong() const { return rand() % 3 - 1; }

  void initial();

  void initial_candidate_loc();
  /**
   * @brief forbid all virtual links
   *
   */
  void forbidVirLinks();
  /**
   * @brief set weight of virtual link to 0
   *
   */
  void releaseVirLinks();

  void lower_update();

  // void supper_update();

  /**
   * @brief when the number of user is small then there are some condition to
   * directly get optimization solution
   *
   *
   * @return ture if find optimization solution, false otherwise
   */
  bool smallestUer();

  bool bestLayoutFlow(Server_loc &server);

  bool more_check(Server_loc &server);

  bool domain_intersection_check();

  /**
 * choose some vertices as servers then compute min cost max flow, if there sum
 * of
 * pre link value of node v greater or equal than server_price then add this v
 * as a new server location
 *
 * @param server
 * @param recursive is recursive call
 */
  void update(Server_loc &server, bool recursive = true);

  void initial_case(void);

  void delete_candiate(void);

  void generate_case(Server_loc &lhs, Server_loc &rhs);

  void randdom_generate(void);

  char *output();

  double start_time;
  Server_loc best_loc;


 public:
  Loc_choose(undirected_graph &g, int network_n_num, int user_n_num,
             int serice_p, int vir_source, int vir_target, int target_vir_links,
             const vector<int> &node_map, const vector<int> &dcaps,
             const vector<int> &caps)
      : graph(g),
        network_node_num(network_n_num),
        user_node_num(user_n_num),
        server_price(serice_p),
        virtual_source(vir_source),
        virtual_target(vir_target),
        target_link_start(target_vir_links),
        network_node_user_map(node_map),
        user_demand(dcaps),
        orignal_caps(caps) {
    start_time = systemTime();
    best_loc.success_bw=-100;
    best_loc.total_price=INF;
    totCap = 0;

    for (vector<int>::const_iterator it = user_demand.begin();
         it != user_demand.end(); it++) {
      totCap += *it;
    }
    
    choose_time.resize(network_node_num, 0  );
    
    server_candiate_locs.resize(user_node_num);
    user_to_network_map.resize(user_node_num);

    value_supper = user_demand.size() * server_price;

    value_lower = INF;
    user_direct_server.resize(user_node_num, false);
    network_to_user_inv_distance.resize(network_node_num);
    sum_of_pass_flow.resize(network_node_num);
    for (int i = 0; i < network_node_num; i++) {
      network_to_user_inv_distance[i].resize(user_node_num, 0.0f);
      sum_of_pass_flow[ i ].resize( network_node_num, 0 );
    }

    if (network_node_num < 300) {
      para.first_class_candiate_num = 4;
    }
    initial();
  }
  char *solve();
};
}

#endif
