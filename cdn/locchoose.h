#ifndef _LOC_CHOOSE_H
#define _LOC_CHOOSE_H

#include<vector>
#include <utility>
#include <limits>
#include<set>
#include<map>
#include "graph.h"

namespace raptor {

using namespace std;




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
  vector<Server_loc> server_candiate;

  int lest_link_price;
  int large_link_price;
  double mean_link_price;
  int middle_link_price;

  int smallest_user;
  int large_user;

  int  randTryNum;
  /**
   * for every user u there at less a network node in server_candiate_locs[u] used as a server
   * 
   */

  vector<vector<int>> server_candiate_locs;

  vector<int> user_to_network_map;

  vector<bool> user_direct_server;

  set<int> choosedServer;

  vector<int> allChoose;
  
  double delete_para;

  int saoDong(  ) const{
    return  rand() % 5 -2;
  }
  
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

  void tryKServer(const int k);

  void bestLayoutFlow(Server_loc &server);

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


  void  delete_canduate( void );

  char *output();
  


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

    randTryNum=100;
    delete_para = 0.5;
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

    initial();
  }
  char *solve();
};
}

#endif
