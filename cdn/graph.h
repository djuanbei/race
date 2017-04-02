#ifndef __graph_h
#define __graph_h

#include <limits>
#include "heap.h"

#include <vector>

namespace raptor {

using namespace std;
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
  VI dist, width, flow, caps, backDis;
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

  void dijkstra_tree(const int src, vector<int> &dis);
    
  void dijkstra_reverse_tree(const int snk, const vector<int> &caps, vector<int> &dis);

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
                                   const vector<int> &caps, const int totCap, vector<int> &outs,
                                   vector<int> &ins,
                                   vector<int> &node_sum_value , vector<vector<int> > &sum_of_pass_flow);

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
}

#endif
