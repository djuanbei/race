#ifndef __ROUTE_H__
#define __ROUTE_H__

#include "lib_io.h"

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <vector>

#include <utility>
namespace raptor {

using namespace std;

class Fixed_heap {
 private:
  vector<pair<int, int>> heap;  // Fixed_heap of Keys
  vector<int> indices;          // Each Key's position (index) in the Fixed_heap

  int size;

  // Index "traversal" functions
  static inline int left(int i) { return (i << 1) | 1; }
  static inline int right(int i) { return (i + 1) << 1; }
  static inline int parent(int i) { return (i - 1) >> 1; }
  void percolateUp(int i);

  void percolateDown(int i);

 public:
  Fixed_heap(int cap) : heap(cap), indices(cap, -1), size(0) {}

  bool empty() const { return 0 == size; }

  void push(pair<int, int> k);

  pair<int, int> &top() { return heap[0]; }

  const pair<int, int> &operator[](const int i) const { return heap[i]; }

  void pop();

  int len() const { return size; }

  void clear() {
    size = 0;
    fill(indices.begin(), indices.end(), -1);
  }
};

const static int INF = numeric_limits<int>::max() / 4;
class undirected_graph {
  typedef vector<int> VI;

  typedef vector<VI> VVI;

  typedef pair<int, int> PII;
  typedef vector<PII> VPII;

  typedef undirected_graph graph_t;
  struct endElement {
    int weight;
    int snk;
    endElement() : weight(0), snk(0) {}
  };

  struct startElement {
    int link;
    int src;
    startElement() : link(0), src(0) {}
  };

  struct tempElment {
    int id;
    int src;
    bool operator<(const tempElment &other) const { return src < other.src; }
  };

 private:
  int vertex_num;
  int link_num;

  vector<int> _srcs;
  // out links
  vector<int> outIndex;
  vector<endElement> link_ends;
  // in links
  vector<int> inIndex;
  vector<startElement> link_starts;

  vector<int> in2outLink_map;  // link_ends  --> orignal link
  vector<int> out2inLink_map;  // orginal link --> link_end
  vector<int> inPeerLinkMap;   // link_end peer link

  vector<bool> found;
  VI dist, width, flow;
  Fixed_heap Q;
  VPII dad;

  int getPeerLink(const int id) const { return id ^ 1; }
  bool _findSrc(const int link, int &src) const {
    src = vertex_num + 1;
    if (link < 0 || link >= link_num) return false;
    src = _srcs[link];
    return true;
  }
  inline bool _findSnk(const int link, int &snk) const {
    snk = vertex_num + 1;
    if (link >= link_num) return false;
    snk = link_ends[link].snk;
    return true;
  }

  inline bool _findSrcSnk(const int link, int &src, int &snk) const {
    if (!_findSrc(link, src)) return false;
    _findSnk(link, snk);
    return true;
  }

  int dijkstra(const int s, const int t, const vector<int> &caps,
               vector<int> &dist, vector<int> &width, vector<int> &flow,
               Fixed_heap &Q);

  int dijkstra(const int s, const int t, const vector<int> &caps);

  int dijkstra(const int src, const int snk, const vector<int> &caps,
               vector<int> &link_path);

 public:
  undirected_graph() : vertex_num(0), link_num(0), Q(10) {}
  void clear() {
    vertex_num = link_num = 0;
    _srcs.clear();
    outIndex.clear();
    link_ends.clear();
    inIndex.clear();
    link_starts.clear();
    in2outLink_map.clear();
    out2inLink_map.clear();
  }

  /**
   * @brief initial undirected graph by undieced links
   *
   * @param esrcs  link source id
   * @param esnks  link target id
   * @param eweights  link weight
   */

  void initial(const vector<int> &esrcs, const vector<int> &esnks,
               const vector<int> &eweights);

  inline size_t getVertex_num(void) const { return vertex_num; }

  inline int getLink_num(void) const { return link_num; }

  inline int getOutDegree(int vertex) const {
    if (vertex >= vertex_num) return 0;
    return outIndex[vertex + 1] - outIndex[vertex];
  }

  inline int getInDegree(int vertex) const {
    if (vertex >= vertex_num) return 0;
    return inIndex[vertex + 1] - inIndex[vertex];
  }

  inline bool findSrc(const int link, int &src) const {
    return _findSrc(out2inLink_map[link], src);
  }
  inline bool findSnk(const int link, int &snk) const {
    return _findSnk(out2inLink_map[link], snk);
  }
  bool findRhs(const int link, const int lhs, int &rhs) const;

  inline bool findSrcSnk(const int link, int &src, int &snk) const {
    return _findSrcSnk(out2inLink_map[link], src, snk);
  }

  inline int getWeight(const int link) const {
    return link_ends[out2inLink_map[link]].weight;
  }

  inline void setWeight(const int link, int w) {
    link_ends[out2inLink_map[link]].weight = w;
  }

  int getAdj(int v, int i) const { return in2outLink_map[outIndex[v] + i]; }
  int getReAdj(int v, int i) const {
    return in2outLink_map[link_starts[inIndex[v] + i].link];
  }

  /**
   *
   * compute shortest path tree which distance less or equal than limit
   * @param src  source vertex id
   * @param limit  distance limit
   * @param tree  pair<int, int> is  pair of vertex id and distance.
   */
  void dijkstra_limit_tree(const int src, const int limit,
                           vector<pair<int, int>> &tree);

  /**
   * greed method to compute approximation minimum vertex cover
   *
   * @param nodes choosed cover set
   */
  void getMinVertexCover(vector<int> &nodes);

  /**
   *@brief compute shortest path which connect src and snk
   *
   * @param src source vertex id
   * @param snk target vertex id
   * @param path __return__ link path
   *
   * @return true if find a path, false otherwise
   */

  bool compute_shortest_path_dijkstra(const int src, const int snk,

                                      vector<int> &path);

  /**
   *@brief compute a shortest path which connect src and snk by bidirection
   *dijstra method
   *
   * @param src source vertex id
   * @param snk target vertex id
   * @param path __return__ link path
   *
   * @return true if find a path, false otherwise
   */

  bool bicompute_shortest_path_dijkstra(const int src, const int snk,
                                        vector<int> &path);

  /**
   *  sequence call
   * compute the mincost maxflow from src to snk
   * @param src
   * @param snk
   * @param caps  undirected link capacity
   * @param outs sum of output bandwidth from every server
   * @param ins sum of every user in bandwidth
   * @param node_sum_value corresponding sum of node value for this min cost max
   * *
   * @return  maxflow  and total cost
   */
  pair<int, int> getMinCostMaxFlow(const int src, const int snk,
                                   const vector<int> &caps, vector<int> &outs,
                                   vector<int> &ins,
                                   vector<int> &node_sum_value);

  /**
   *  parallel call
   * compute the mincost maxflow from src to snk
   * @param src
   * @param snk
   * @param caps undirected link capacity
   *
   * @return maxflow  and total cost
   */

  pair<int, int> getMinCostMaxFlowP(const int src, const int snk,
                                    const vector<int> &caps);

  /**
   *
   *
   * @param src  the source vertex id
   * @param snk  the target vertex id
   * @param caps   undirected link capacity
   * @param node_paths  __return__ min cost malfow split node paths
   * @param bws  the bandwidth of every split flow

   * flow
   */
  void getMinCostMaxFlow(int src, const int snk, const vector<int> &caps,
                         vector<vector<int>> &node_paths, vector<int> &bws);

  /**
   * @brief check wether this path is connect src and snk
   *
   * @param src  source vertex id
   * @param snk  target vertex id
   * @param path  link path
   *
   * @return  true if path connects src and snk, false otherwise
   */
  bool isValidatePath(const int &src, const int &snk,
                      const vector<int> &path) const;
  /**
   *
   *
   * @param path  link path
   *
   * @return  sum of weigth in this path
   */
  int path_cost(const vector<int> &path) const;
};

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

  /**
   * for every user u there at less a network node in server_candiate_locs[u] used as a server
   * 
   */

  vector<vector<int>> server_candiate_locs;

  vector<int> user_to_network_map;

  vector<bool> user_direct_server;

  set<int> choosedServer;

  vector<int> allChoose;

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

  

  char *output();
  double delete_para;

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
    delete_para = 0.5;
    totCap = 0;

    for (vector<int>::const_iterator it = user_demand.begin();
         it != user_demand.end(); it++) {
      totCap += *it;
    }
    server_candiate_locs.resize(user_node_num);
    user_to_network_map.resize(user_node_num);

    value_supper = user_demand.size() * server_price;

    value_lower = 0;
    user_direct_server.resize(user_node_num, false);

    initial();
  }
  char *solve();
};
}
void deploy_server(char *graph[MAX_EDGE_NUM], int edge_num, char *filename);

#endif
