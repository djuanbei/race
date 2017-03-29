#ifndef _LOC_CHOOSE_H
#define _LOC_CHOOSE_H

#include <cstdlib>
#include <limits>
#include <map>
#include <set>
#include <utility>
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

    Server_loc() : success_bw(0), total_price(0) {}
    bool operator<(const Server_loc &other) const {
      
      // if (success_bw > other.success_bw) {
      //   return true;
      // } else if (success_bw < other.success_bw) {
      //   return false;
      // }
      

      if (0 == success_bw) {
        return 0 == other.success_bw;
      }
      if (0 == other.success_bw) {
        return true;
      }
      return (total_price / (success_bw + 0.01)) <
             (other.total_price / (other.success_bw + 0.01));
    }
  };



  struct Para{
    
    int randTryNum;
    
    double delete_para;

    int first_class_candiate_num;
    
    bool deleteSmall;
    
    int stable_bound;
    int success_left_num;
    int part_left_num;
    int add_num;
    
    
    Para(  ){
      
      randTryNum= 100;
      delete_para = 0.5;
      first_class_candiate_num = 1;
      deleteSmall = false;
      stable_bound=26;
      success_left_num=10;
      part_left_num=50;
      add_num=7;
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
  
  vector<Server_loc> part_server_candiate;

  int lest_link_price;
  
  int large_link_price;
  double mean_link_price;
  int middle_link_price;

  int smallest_user;
  int large_user;


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
  
  vector<vector<float>> network_to_user_inv_distance;
  
  Para para;

  void addLoc( Server_loc & server ){
    if( totCap==server.success_bw ){
      success_server_candiate.push_back( server );
    }
    part_server_candiate.push_back( server );

  }
  
  bool time_end(  ) const{
    return systemTime() - start_time > time_bound - 10;
  }
  
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

  void supper_update();

  /**
   * @brief when the number of user is small then there are some condition to
   * directly get optimization solution
   *
   *
   * @return ture if find optimization solution, false otherwise
   */
  bool smallestUer();


  bool bestLayoutFlow(Server_loc &server);

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

  void delete_canduate(void);

  void generate_case(Server_loc &lhs, Server_loc &rhs);

  void randdom_generate(void);

  char *output();

  double start_time;

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
    
    totCap = 0;

    for (vector<int>::const_iterator it = user_demand.begin();
         it != user_demand.end(); it++) {
      totCap += *it;
    }
    server_candiate_locs.resize(user_node_num);
    user_to_network_map.resize(user_node_num);

    value_supper = user_demand.size() * server_price;

    value_lower = INF;
    user_direct_server.resize(user_node_num, false);
    network_to_user_inv_distance.resize(network_node_num);
    for (int i = 0; i < network_node_num; i++) {
      network_to_user_inv_distance[i].resize(user_node_num, 0.0f);
    }

    initial();
  }
  char *solve();
};
}

#endif
