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

namespace raptor {

using namespace std;

class Fixed_heap {
private:
  vector<pair<int, int>> heap; // Fixed_heap of Keys
  vector<int> indices;         // Each Key's position (index) in the Fixed_heap

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

class undirected_graph {

  typedef vector<int> VI;

  typedef vector<VI> VVI;

  typedef pair<int, int> PII;
  typedef vector<PII> VPII;

  const int INF = numeric_limits<int>::max() / 4;

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

  vector<int> in2outLink_map; // link_ends  --> orignal link
  vector<int> out2inLink_map; // orginal link --> link_end

  vector<bool> found;
  VI dist, pi, width, flow;
  Fixed_heap Q;
  VPII dad;

  int getPeerLink(const int id) const { return id ^ 1; }
  bool _findSrc(const int link, int &src) const {
    src = vertex_num + 1;
    if (link < 0 || link >= link_num)
      return false;
    src = _srcs[link];
    return true;
  }
  inline bool _findSnk(const int link, int &snk) const {
    snk = vertex_num + 1;
    if (link >= link_num)
      return false;
    snk = link_ends[link].snk;
    return true;
  }

  inline bool _findSrcSnk(const int link, int &src, int &snk) const {
    if (!_findSrc(link, src))
      return false;
    _findSnk(link, snk);
    return true;
  }

  int dijkstra(const int s, const int t, const vector<int> &caps,
               vector<int> &dist, vector<int> &width, vector<int> &flow,
               vector<int> &pi, Fixed_heap &Q);

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

  void initial(const vector<int> &esrcs, const vector<int> &esnks,
               const vector<int> &eweights);

  inline size_t getVertex_num(void) const { return vertex_num; }

  inline int getLink_num(void) const { return link_num; }

  inline int getOutDegree(int vertex) const {
    if (vertex >= vertex_num)
      return 0;
    return outIndex[vertex + 1] - outIndex[vertex];
  }

  inline int getInDegree(int vertex) const {
    if (vertex >= vertex_num)
      return 0;
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

  /**
   *
   *
   * @param src
   * @param snk
   * @param path
   *
   * @return true if find a path, false otherwise
   */

  bool compute_shortest_path_dijkstra(const int src, const int snk,

                                      vector<int> &path);

  /**
   *
   *
   * @param src
   * @param snk
   * @param path
   *
   * @return true if find a path, false otherwise
   */

  bool bicompute_shortest_path_dijkstra(const int src, const int snk,
                                        vector<int> &path);

  /** 
   *  sequence call
   * 
   * @param src 
   * @param snk 
   * @param caps  undirected link capacity
   * 
   * @return  maxflow  and total cost
   */
  pair<int, int> getMinCostMaxFlow(const int src, const int snk,
                                   const vector<int> &caps);

  /** 
   *  parallel call
   * 
   * @param src 
   * @param snk 
   * @param caps undirected link capacity
   * 
   * @return maxflow  and total cost
   */
  pair<int, int> getMinCostMaxFlowP(const int src, const int snk,
                                    const vector<int> &caps);
  void getMinCostMaxFlow(int src, const int snk, const vector<int> &caps,
                         vector<vector<int>> &node_paths, vector<int> &bws);

  bool isValidatePath(const int &src, const int &snk,
                      const vector<int> &path) const;

  int path_cost(const vector<int> &path) const;
};
}
void deploy_server(char *graph[MAX_EDGE_NUM], int edge_num, char *filename);

#endif
