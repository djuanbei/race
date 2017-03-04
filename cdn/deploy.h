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


using namespace std;


class Fixed_heap {
 private:
  vector<pair<int, int>> heap;  // Fixed_heap of Keys
  vector<int> indices;      // Each Key's position (index) in the Fixed_heap

  int size;

  // Index "traversal" functions
  static inline int left(int i) { return (i << 1) | 1; }
  static inline int right(int i) { return (i + 1) << 1; }
  static inline int parent(int i) { return (i - 1) >> 1; }
  void percolateUp(int i) {
    pair<int, int> x = heap[i];
    int p = parent(i);

    while (i != 0 && (x.first< heap[p].first)) {
      heap[i] = heap[p];
      indices[heap[i].second] = i;
      i = p;
      p = parent(p);
    }
    heap[i] = x;
    indices[x.second] = i;
  }

  void percolateDown(int i) {
    pair<int, int> x = heap[i];
    while (left(i) < size) {
      int child = right(i) < size && (heap[right(i)].first< heap[left(i)].first)
                      ? right(i)
                      : left(i);
      if (!(heap[child].first< x.first)) break;
      heap[i] = heap[child];
      indices[heap[i].second] = i;
      i = child;
    }
    heap[i] = x;
    indices[x.second] = i;
  }

 public:
  Fixed_heap( int cap)
      : heap(cap), indices(cap, -1),  size(0) {}

  bool empty() const { return 0 == size; }
  void push(pair<int, int> k) {
    if (-1 == indices[k.second]) {
      heap[size] = k;
      indices[k.second] = size;
      percolateUp(size);
      size++;
    } else {
      heap[indices[k.second]].first = k.first;
      percolateUp(indices[k.second]);
    }
  }
  pair<int, int> &top() { return heap[0]; }

  const pair<int, int> &operator[](const int i) const { return heap[i]; }

  void pop() {
    indices[heap[0].second] = -1;
    heap[0] = heap[size - 1];
    if (size > 0) {
      indices[heap[0].second] = 0;
      size--;
      percolateDown(0);
    }
  }

  int len() const { return size; }

  void clear() {
    size = 0;
    fill(indices.begin(), indices.end(), -1);
  }
};



class directed_graph {
  typedef pair<int, int> PII;
  typedef  directed_graph graph_t;
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
    bool operator <(const tempElment &other ) const{
      return src< other.src;
    }
  };


 private:
  bool save_path;
  int infi_value;

  size_t ksp_k;

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


 public:

  directed_graph(const size_t k = 1)
      : save_path(true),
        infi_value(numeric_limits<int>::max() / 10e10),
        ksp_k(k),
        vertex_num(0),
        link_num(0) {
    assert(k > 0);
  }
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

  void initial(vector<int> &srcs, vector<int> &snks, vector<int> &weights,
               bool save = false) {
    assert(srcs.size() == snks.size() && srcs.size() == weights.size());
    clear();
    save_path = save;

    if (0 == srcs.size()) return;
    link_num = srcs.size();
    in2outLink_map.resize(link_num);
    out2inLink_map.resize(link_num);

    vector<tempElment> tContian;

    int i, j, ver_max;
    tempElment dummy_pred;

    ver_max = srcs[0];

    for (i = 0; i < link_num; i++) {
      if (ver_max < srcs[i]) {
        ver_max = srcs[i];
      }
      if (ver_max < snks[i]) {
        ver_max = snks[i];
      }
    }
    vertex_num = ver_max + 1;
    if (vertex_num > 5000) save_path = false;



    for (size_t i = 0; i < srcs.size(); i++) {
      dummy_pred.id = i;
      dummy_pred.src = srcs[i];
      tContian.push_back(dummy_pred);
    }
    std::sort(tContian.begin(), tContian.end());

    outIndex.push_back(0);
    i = 0;
    while (i < tContian[0].src) {
      outIndex.push_back(0);
      i++;
    }

    endElement temp;
    temp.weight = weights[tContian[0].id];
    temp.snk = snks[tContian[0].id];

    _srcs.push_back(srcs[tContian[0].id]);
    link_ends.push_back(temp);
    in2outLink_map[0] = tContian[0].id;
    out2inLink_map[tContian[0].id] = 0;
    i = 1;
    for (; i < (int)tContian.size(); i++) {
      if (tContian[i].src != tContian[i - 1].src) {
        for (j = tContian[i - 1].src; j < tContian[i].src; j++) {
          outIndex.push_back(link_ends.size());
        }
      }
      temp.weight = weights[tContian[i].id];
      temp.snk = snks[tContian[i].id];

      _srcs.push_back(srcs[tContian[i].id]);
      link_ends.push_back(temp);
      in2outLink_map[i] = tContian[i].id;
      out2inLink_map[tContian[i].id] = i;
    }

    for (j = tContian[i - 1].src; j < vertex_num; j++) {
      outIndex.push_back(link_ends.size());
    }

    tContian.clear();
    for (size_t i = 0; i < srcs.size(); i++) {
      dummy_pred.id = i;
      dummy_pred.src = snks[i];
      tContian.push_back(dummy_pred);
    }

    std::sort(tContian.begin(), tContian.end());

    inIndex.push_back(0);
    i = 0;
    while (i < tContian[0].src) {
      inIndex.push_back(0);
      i++;
    }
    startElement dummy;

    dummy.link = out2inLink_map[tContian[0].id];
    dummy.src = srcs[tContian[0].id];

    link_starts.push_back(dummy);
    i = 1;
    for (; i < (int)tContian.size(); i++) {
      if (tContian[i].src != tContian[i - 1].src) {
        for (j = tContian[i - 1].src; j < tContian[i].src; j++) {
          inIndex.push_back(link_starts.size());
        }
      }

      dummy.link = out2inLink_map[tContian[i].id];
      dummy.src = srcs[tContian[i].id];
      link_starts.push_back(dummy);
    }

    for (j = tContian[i - 1].src; j < vertex_num; j++) {
      inIndex.push_back(link_starts.size());
    }

    for (size_t i = 0; i < srcs.size(); i++) {
      int temp;
      findSrc(i, temp);
      assert(temp == srcs[i]);
    }
  }

  bool isDirect() const { return true; }

  inline size_t getVertex_num(void) const { return vertex_num; }

  void setInfi(const int infi) { infi_value = infi; }

  inline int getLink_num(void) const { return link_num; }

  inline int getOutDegree(int vertex) const {
    if (vertex >= vertex_num) return 0;
    return outIndex[vertex + 1] - outIndex[vertex];
  }

  inline int getInDegree(int vertex) const {
    if (vertex >= vertex_num) return 0;
    return inIndex[vertex + 1] - inIndex[vertex];
  }

  inline endElement &getLink(int vertex, int k) {
    assert(vertex < vertex_num && k < getOutDegree(vertex));
    return link_ends[outIndex[vertex] + k];
  }

  inline int getWeight(const int link) const {
    return link_ends[out2inLink_map[link]].weight;
  }

  int getAdj(int v, int i) const { return in2outLink_map[outIndex[v] + i]; }

  int getReAdj(int v, int i) const {
    return in2outLink_map[link_starts[inIndex[v] + i].link];
  }

  inline bool findSrc(const int link, int &src) const {
    return _findSrc(out2inLink_map[link], src);
  }
  inline bool findSnk(const int link, int &snk) const {
    return _findSnk(out2inLink_map[link], snk);
  }
  bool findRhs(const int link, const int lhs, int &rhs) const {
    int tempSrc, tempSnk;
    if (!findSrcSnk(link, tempSrc, tempSnk)) {
      return false;
    }
    if (lhs == tempSrc) {
      rhs = tempSnk;
      return true;
    }
    if (lhs == tempSnk) {
      rhs = tempSrc;
      return true;
    }

    return false;
  }

  inline bool findSrcSnk(const int link, int &src, int &snk) const {
    return _findSrcSnk(out2inLink_map[link], src, snk);
  }

  inline void setLinkWeight(const int link, const int weight) {
    link_ends[out2inLink_map[link]].weight = weight;
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

                                      vector<int> &path) {
    path.clear();
    if (src >= vertex_num || snk >= vertex_num) return false;
    if (src == snk) return true;
    size_t j, outDegree, link, next;
    int current;
    int weight;
    vector<int> dis(vertex_num, infi_value);
    vector<int> preLink(vertex_num, link_num + 1);
    vector<int> parent(vertex_num, -1);
    dis[src] = 0;
    Fixed_heap Q( vertex_num);

    Q.push(make_pair(0.0, src));

    while (!Q.empty()) {
      PII p = Q.top();
      current = p.second;
      if (current == snk) {
        while (current != src) {
          path.push_back(in2outLink_map[preLink[current]]);

          current = _srcs[preLink[current]];
        }
        reverse(path.begin(), path.end());
        return true;
      }
      Q.pop();

      current = p.second;
      outDegree = getOutDegree(current);
      int current_weight = p.first;
      for (j = 0; j < outDegree; j++) {
        link = outIndex[current] + j;
        const endElement &neighbour = link_ends[link];
        weight = current_weight + neighbour.weight;
        next = neighbour.snk;
        if (weight < dis[snk] && weight < dis[next]) {
          parent[next] = current;
          preLink[next] = link;
          dis[next] = weight;
          Q.push(make_pair(weight, next));
        }
      }
    }
    return false;
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

  bool bicompute_shortest_path_dijkstra(const int src, const int snk,
                                        vector<int> &path) {
    path.clear();
    if (src >= vertex_num || snk >= vertex_num) return false;
    if (src == snk) return true;

    size_t j, outDegree, inDegree, link, next;
    int current;
    int weight;

    vector<int> dis(vertex_num, infi_value);
    vector<int> preLink(vertex_num, -1);
    vector<bool> check(vertex_num, false);
    dis[src] = 0;

    Fixed_heap Q( vertex_num);

    Q.push(make_pair(0.0, src));

    vector<int> bdis(vertex_num, infi_value);
    vector<int> bpreLink(vertex_num, -1);
    bdis[snk] = 0;
    Fixed_heap bQ( vertex_num);

    bQ.push(make_pair(0.0, snk));

    bool re = false;
    int best_dis = infi_value;
    int best_current = -1;

    while (!Q.empty() && !bQ.empty()) {
      PII p = Q.top();
      current = p.second;
      if (check[current]) {
        re = true;
        break;
      }
      check[current] = true;
      Q.pop();
      outDegree = getOutDegree(current);
      int current_weight = dis[current];
      for (j = 0; j < outDegree; j++) {
        link = outIndex[current] + j;
        const endElement &neighbour = link_ends[link];
        weight = current_weight + neighbour.weight;
        next = neighbour.snk;

        if (weight < dis[next]) {
          preLink[next] = link;
          dis[next] = weight;
          Q.push(make_pair(weight, next));
        }
      }

      p = bQ.top();
      current = p.second;
      if (check[current]) {
        re = true;
        break;
      }

      check[current] = true;
      bQ.pop();
      inDegree = getInDegree(current);
      current_weight = bdis[current];
      for (j = 0; j < inDegree; j++) {
        link = inIndex[current] + j;

        const startElement &neighbour = link_starts[link];
        link = neighbour.link;
        weight = current_weight + link_ends[link].weight;
        next = neighbour.src;

        if (weight < bdis[next]) {
          bpreLink[next] = link;
          bdis[next] = weight;
          bQ.push(make_pair(weight, next));
        }
      }
    }

    if (re) {
      int temp_dis;
      int temp;
      int num = Q.len();
      for (int i = 0; i < num; i++) {
        temp = Q[i].second;
        temp_dis = dis[temp] + bdis[temp];
        if (temp_dis < best_dis) {
          best_dis = temp_dis;
          best_current = temp;
        }
      }
      num = bQ.len();
      for (int i = 0; i < num; i++) {
        temp = bQ[i].second;
        temp_dis = dis[temp] + bdis[temp];
        if (temp_dis < best_dis) {
          best_dis = temp_dis;
          best_current = temp;
        }
      }

      current = best_current;

      while (current != src) {
        path.push_back(in2outLink_map[preLink[current]]);
        current = _srcs[preLink[current]];
      }
      reverse(path.begin(), path.end());

      temp = best_current;
      while (temp != snk) {
        path.push_back(in2outLink_map[bpreLink[temp]]);
        _findSnk(bpreLink[temp], temp);
      }
    }

    return re;
  }



  bool isValidatePath(const int &src, const int &snk, const vector<int> &path) const {
    int current = src;
    int next = src;
    int next1;
    for (typename vector<int>::const_iterator it = path.begin(); it != path.end();
         it++) {
      if (!findSrcSnk(*it, current, next1)) {
        return false;
      }

      if (current != next) {
        return false;
      }
      next = next1;
    }

    return next == snk;
  }
  int path_cost(const vector<int> &path) const {
    int re = 0;
    for (typename vector<int>::const_iterator it = path.begin(); it != path.end();
         it++) {
      re += link_ends[out2inLink_map[*it]].weight;
    }
    return re;
  }

  static void printPath(const vector<int> &path) {
    if (0 == path.size()) return;

    size_t i = 0;
    string ss;

    while (i < path.size()) {
      std::cout << "by " << path[i] << std::endl;
      i++;
    }
  }
};


class undirected_graph {
  struct startElement {
    int link;
    int src;
    startElement() : link(0), src(0) {}
  };
  struct tempElment {
    int id;
    int src;
    bool operator < ( const tempElment &other  ) const{
      return src< other.src;
    }
  };



 private:
  int vertex_num;
  int link_num;

  vector<int> _srcs;
  // out links
  vector<int> outIndex;
  vector<int> link_ends;
  // in links
  vector<int> inIndex;
  vector<startElement> link_starts;

  vector<int> in2outLink_map;  // link_ends  --> orignal link
  vector<int> out2inLink_map;  // orginal link --> link_end

  bool _findSrc(const int link, int &src) const {
    src = vertex_num + 1;
    if (link < 0 || link >= link_num) return false;
    src = _srcs[link];
    return true;
  }
  inline bool _findSnk(const int link, int &snk) const {
    snk = vertex_num + 1;
    if (link >= link_num) return false;
    snk = link_ends[link];
    return true;
  }
  inline bool _findSrcSnk(const int link, int &src, int &snk) const {
    if (!_findSrc(link, src)) return false;
    _findSnk(link, snk);
    return true;
  }

 public:
  void initial(vector<int> &srcs, vector<int> &snks) {
    if (0 == srcs.size()) return;
    link_num = srcs.size();
    in2outLink_map.resize(link_num);
    out2inLink_map.resize(link_num);

    vector<tempElment> tContian;

    int i, j, ver_max;
    tempElment dummy_pred;

    ver_max = srcs[0];

    for (i = 0; i < link_num; i++) {
      if (ver_max < srcs[i]) {
        ver_max = srcs[i];
      }
      if (ver_max < snks[i]) {
        ver_max = snks[i];
      }
    }
    vertex_num = ver_max + 1;

    for (size_t i = 0; i < srcs.size(); i++) {
      dummy_pred.id = i;
      dummy_pred.src = srcs[i];
      tContian.push_back(dummy_pred);
    }
    std::sort(tContian.begin(), tContian.end());

    outIndex.push_back(0);
    i = 0;
    while (i < tContian[0].src) {
      outIndex.push_back(0);
      i++;
    }

    int temp;

    temp = snks[tContian[0].id];

    _srcs.push_back(srcs[tContian[0].id]);
    link_ends.push_back(temp);
    in2outLink_map[0] = tContian[0].id;
    out2inLink_map[tContian[0].id] = 0;
    i = 1;
    for (; i < (int)tContian.size(); i++) {
      if (tContian[i].src != tContian[i - 1].src) {
        for (j = tContian[i - 1].src; j < tContian[i].src; j++) {
          outIndex.push_back(link_ends.size());
        }
      }

      temp = snks[tContian[i].id];

      _srcs.push_back(srcs[tContian[i].id]);
      link_ends.push_back(temp);
      in2outLink_map[i] = tContian[i].id;
      out2inLink_map[tContian[i].id] = i;
    }

    for (j = tContian[i - 1].src; j < vertex_num; j++) {
      outIndex.push_back(link_ends.size());
    }

    tContian.clear();
    for (size_t i = 0; i < srcs.size(); i++) {
      dummy_pred.id = i;
      dummy_pred.src = snks[i];
      tContian.push_back(dummy_pred);
    }

    std::sort(tContian.begin(), tContian.end());

    inIndex.push_back(0);
    i = 0;
    while (i < tContian[0].src) {
      inIndex.push_back(0);
      i++;
    }
    startElement dummy;

    dummy.link = out2inLink_map[tContian[0].id];
    dummy.src = srcs[tContian[0].id];

    link_starts.push_back(dummy);
    i = 1;
    for (; i < (int)tContian.size(); i++) {
      if (tContian[i].src != tContian[i - 1].src) {
        for (j = tContian[i - 1].src; j < tContian[i].src; j++) {
          inIndex.push_back(link_starts.size());
        }
      }

      dummy.link = out2inLink_map[tContian[i].id];
      dummy.src = srcs[tContian[i].id];
      link_starts.push_back(dummy);
    }

    for (j = tContian[i - 1].src; j < vertex_num; j++) {
      inIndex.push_back(link_starts.size());
    }
  }

  bool isDirect() const { return false; }
  inline size_t getVertex_num(void) const { return vertex_num; }
  inline int getLink_num(void) const { return link_num; }

  inline int getOutDegree(int vertex) const {
    if (vertex >= vertex_num) return 0;
    return outIndex[vertex + 1] - outIndex[vertex] + inIndex[vertex + 1] -
           inIndex[vertex];
  }

  inline int getInDegree(int vertex) const {
    if (vertex >= vertex_num) return 0;
    return outIndex[vertex + 1] - outIndex[vertex] + inIndex[vertex + 1] -
           inIndex[vertex];
  }
  int getAdj(int vertex, int i) const {
    int k = outIndex[vertex + 1] - outIndex[vertex];
    if (i < k)
      return in2outLink_map[outIndex[vertex] + i];
    else
      return in2outLink_map[link_starts[inIndex[vertex] + i - k].link];
  }
  int getReAdj(int vertex, int i) const {
    int k = inIndex[vertex + 1] - inIndex[vertex];
    if (i < k)
      return in2outLink_map[link_starts[inIndex[vertex] + i].link];
    else
      return in2outLink_map[outIndex[vertex] + i - k];
  }
  bool findRhs(const int link, const int lhs, int &rhs) const {
    int tempSrc, tempSnk;
    if (!findSrcSnk(link, tempSrc, tempSnk)) {
      return false;
    }
    if (lhs == tempSrc) {
      rhs = tempSnk;
      return true;
    }
    if (lhs == tempSnk) {
      rhs = tempSrc;
      return true;
    }

    return false;
  }
  inline bool findSrcSnk(const int link, int &src, int &snk) const {
    return _findSrcSnk(out2inLink_map[link], src, snk);
  }

  bool findSrc(const int link, int &src) const {
    return _findSrc(out2inLink_map[link], src);
  }
  inline bool findSnk(const int link, int &snk) const {
    return _findSnk(out2inLink_map[link], snk);
  }
};


void deploy_server(char * graph[MAX_EDGE_NUM], int edge_num, char * filename);

	

#endif
