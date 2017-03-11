#include "deploy.h"

#include "lib_io.h"
#include "lib_time.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

using namespace std;

namespace raptor {

void Fixed_heap::percolateUp(int i) {
  pair<int, int> x = heap[i];
  int p = parent(i);

  while (i != 0 && (x.first < heap[p].first)) {
    heap[i] = heap[p];
    indices[heap[i].second] = i;
    i = p;
    p = parent(p);
  }
  heap[i] = x;
  indices[x.second] = i;
}

void Fixed_heap::percolateDown(int i) {
  pair<int, int> x = heap[i];
  while (left(i) < size) {
    int child = right(i) < size && (heap[right(i)].first < heap[left(i)].first)
                    ? right(i)
                    : left(i);
    if (!(heap[child].first < x.first))
      break;
    heap[i] = heap[child];
    indices[heap[i].second] = i;
    i = child;
  }
  heap[i] = x;
  indices[x.second] = i;
}

void Fixed_heap::push(pair<int, int> k) {
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

void Fixed_heap::pop() {
  indices[heap[0].second] = -1;
  heap[0] = heap[size - 1];
  if (size > 0) {
    indices[heap[0].second] = 0;
    size--;
    percolateDown(0);
  }
}

void undirected_graph::initial(const vector<int> &esrcs,
                               const vector<int> &esnks,
                               const vector<int> &eweights) {
  assert(esrcs.size() == esnks.size() && esrcs.size() == eweights.size());
  clear();

  if (0 == esrcs.size()) {
    return;
  }

  vector<int> srcs, snks, weights;
  for (size_t i = 0; i < esrcs.size(); i++) {
    srcs.push_back(esrcs[i]);
    snks.push_back(esnks[i]);
    weights.push_back(eweights[i]);

    srcs.push_back(esnks[i]);
    snks.push_back(esrcs[i]);
    weights.push_back(eweights[i]);
  }

  link_num = srcs.size();
  in2outLink_map.resize(link_num);
  out2inLink_map.resize(link_num);
  inPeerLinkMap.resize(link_num);
  flow.resize(link_num, 0);

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

  dist.resize(vertex_num, INF);

  width.resize(vertex_num, 0);
  dad.resize(vertex_num);
  Q = Fixed_heap(vertex_num);

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
    temp = in2outLink_map[i];
    temp = getPeerLink(temp);
    inPeerLinkMap[i] = out2inLink_map[temp];

    findSrc(i, temp);
    assert(temp == srcs[i]);
  }
}

bool undirected_graph::findRhs(const int link, const int lhs, int &rhs) const {
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

void undirected_graph::dijkstra_limit_tree(const int src, const int limit,
                                           vector<pair<int, int>> &tree) {

  tree.clear();

  size_t j, outDegree, link, next;
  int current;
  int weight;
  vector<int> dis(vertex_num, INF);
  vector<int> preLink(vertex_num, link_num + 1);
  vector<int> parent(vertex_num, -1);
  dis[src] = 0;
  Fixed_heap Q(vertex_num);

  Q.push(make_pair(0, src));

  while (!Q.empty()) {
    PII p = Q.top();
    current = p.second;

    Q.pop();

    outDegree = getOutDegree(current);
    int current_weight = p.first;
    for (j = 0; j < outDegree; j++) {
      link = outIndex[current] + j;
      const endElement &neighbour = link_ends[link];
      weight = current_weight + neighbour.weight;
      next = neighbour.snk;
      if (weight <= limit && weight < dis[next]) {
        parent[next] = current;
        preLink[next] = link;
        dis[next] = weight;
        Q.push(make_pair(weight, next));
      }
    }
  }

  for (int i = 0; i < vertex_num; i++) {
    if (dis[i] < INF) {
      tree.push_back(make_pair(i, dis[i]));
    }
  }
}

void undirected_graph::getMinVertexCover(vector<int> &nodes) {
  int current, j, outDegree, next, link;
  nodes.clear();
  Fixed_heap Q(vertex_num);
  vector<int> degrees(vertex_num);
  for (int i = 0; i < vertex_num; i++) {
    int outDegree = getOutDegree(i);
    Q.push(make_pair(-outDegree, i));
    degrees[i] = -outDegree;
  }

  vector<bool> check(link_num, false);

  while (!Q.empty()) {
    PII p = Q.top();
    current = p.second;
    Q.pop();
    check[current] = true;

    if (degrees[current] >= 0) {
      continue;
    }

    nodes.push_back(current);
    outDegree = getOutDegree(current);

    for (j = 0; j < outDegree; j++) {
      link = outIndex[current] + j;
      const endElement &neighbour = link_ends[link];
      next = neighbour.snk;
      if (!check[next]) {
        degrees[next]++;
        Q.push(make_pair(degrees[next], next));
      }
    }
  }
}

bool undirected_graph::compute_shortest_path_dijkstra(const int src,
                                                      const int snk,

                                                      vector<int> &path) {
  path.clear();
  if (src >= vertex_num || snk >= vertex_num) {
    return false;
  }

  if (src == snk) {
    return true;
  }

  size_t j, outDegree, link, next;
  int current;
  int weight;
  vector<int> dis(vertex_num, INF);
  vector<int> preLink(vertex_num, link_num + 1);
  vector<int> parent(vertex_num, -1);
  dis[src] = 0;
  Fixed_heap Q(vertex_num);

  Q.push(make_pair(0, src));

  while (!Q.empty()) {
    PII p = Q.top();
    current = p.second;
    if (current == snk) {
      while (current != src) {
        path.push_back(in2outLink_map[preLink[current]] / 2);

        current = _srcs[preLink[current]];
      }
      reverse(path.begin(), path.end());
      return true;
    }
    Q.pop();

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

bool undirected_graph::bicompute_shortest_path_dijkstra(const int src,
                                                        const int snk,
                                                        vector<int> &path) {
  path.clear();
  if (src >= vertex_num || snk >= vertex_num) {
    return false;
  }

  if (src == snk) {
    return true;
  }

  size_t j, outDegree, inDegree, link, next;
  int current;
  int weight;

  vector<int> dis(vertex_num, INF);
  vector<int> preLink(vertex_num, -1);
  vector<bool> check(vertex_num, false);
  dis[src] = 0;

  Fixed_heap Q(vertex_num);

  Q.push(make_pair(0, src));

  vector<int> bdis(vertex_num, INF);
  vector<int> bpreLink(vertex_num, -1);
  bdis[snk] = 0;
  Fixed_heap bQ(vertex_num);

  bQ.push(make_pair(0, snk));

  bool re = false;
  int best_dis = INF;
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
    for (vector<int>::iterator it = path.begin(); it != path.end(); it++) {
      *it /= 2;
    }
  }

  return re;
}

int undirected_graph::dijkstra(const int src, const int snk,
                               const vector<int> &caps, vector<int> &dist,
                               vector<int> &width, vector<int> &flow,
                               Fixed_heap &Q) {

  int j, outDegree, link, next;
  int current;
  int weight;

  fill(dist.begin(), dist.end(), INF);
  fill(width.begin(), width.end(), 0);
  dist[src] = 0;
  width[src] = INF;

  Q.clear();
  Q.push(make_pair(0, src));

  while (!Q.empty()) {
    PII p = Q.top();
    current = p.second;
    if (current == snk) {

      return width[snk];
    }
    Q.pop();

    outDegree = getOutDegree(current);

    int current_weight = p.first;

    for (int j = 0; j < outDegree; j++) {
      link = outIndex[current] + j;
      int cap = (caps[link] - flow[link]);
      if (cap > 0) {
        const endElement &neighbour = link_ends[link];

        next = neighbour.snk;

        weight = dist[current] + neighbour.weight;
        if (weight < dist[snk] && weight < dist[next]) {

          dist[next] = weight;
          dad[next] = make_pair(current, 2 * link + 1);
          width[next] = min(cap, width[current]);
          Q.push(make_pair(weight, next));
        }
      }

      link = inPeerLinkMap[link];

      cap = (flow[link]);
      if (cap > 0) {
        const endElement &neighbour = link_ends[link];

        next = _srcs[link];
        weight = dist[current] - neighbour.weight;

        if (weight < dist[snk] && weight < dist[next]) {

          dist[next] = weight;
          dad[next] = make_pair(current, 2 * link);
          width[next] = min(cap, width[current]);
          Q.push(make_pair(weight, next));
        }
      }
    }
  }

  return 0;
}
int undirected_graph::dijkstra(const int src, const int snk,
                               const vector<int> &caps) {
  return dijkstra(src, snk, caps, dist, width, flow, Q);
}

int undirected_graph::dijkstra(const int src, const int snk,
                               const vector<int> &caps,
                               vector<int> &link_path) {

  link_path.clear();

  size_t j, outDegree, link, next;
  int current;
  int weight;
  vector<int> dis(vertex_num, INF);
  vector<int> preLink(vertex_num, link_num + 1);
  vector<int> parent(vertex_num, -1);
  dis[src] = 0;
  Fixed_heap Q(vertex_num);

  Q.push(make_pair(0, src));

  while (!Q.empty()) {

    PII p = Q.top();
    current = p.second;
    if (current == snk) {
      while (current != src) {
        link_path.push_back(preLink[current]);

        current = _srcs[preLink[current]];
      }
      reverse(link_path.begin(), link_path.end());
      int re = INF;
      for (vector<int>::const_iterator it = link_path.begin();
           it != link_path.end(); it++) {
        re = min(re, caps[*it]);
      }

      return re;
    }
    Q.pop();

    outDegree = getOutDegree(current);
    int current_weight = p.first;
    for (j = 0; j < outDegree; j++) {
      link = outIndex[current] + j;
      if (caps[link] > 0) {
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
  }

  return 0;
}

pair<int, int> undirected_graph::getMinCostMaxFlow(const int src, const int snk,
                                                   const vector<int> &ecaps) {
  assert(link_num == 2 * ecaps.size());
  fill(dist.end(), dist.end(), INF);
  fill(width.begin(), width.end(), 0);
  fill(flow.begin(), flow.end(), 0);

  vector<int> caps(link_num);
  for (size_t i = 0; i < ecaps.size(); i++) {
    int link = out2inLink_map[2 * i];
    caps[link] = ecaps[i];
    link = out2inLink_map[2 * i + 1];
    caps[link] = ecaps[i];
  }

  int totflow = 0, totcost = 0;
  int amt = dijkstra(src, snk, caps);
  while (amt > 0) {
    totflow += amt;
    for (int x = snk; x != src; x = dad[x].first) {
      int link = (dad[x].second) / 2;

      if (1 == ((dad[x].second) % 2)) {
        flow[link] += amt;
        totcost += amt * link_ends[link].weight;
      } else {
        flow[link] -= amt;
        totcost -= amt * link_ends[link].weight;
      }
    }
    amt = dijkstra(src, snk, caps);
  }

  return make_pair(totflow, totcost);
}

pair<int, int> undirected_graph::getMinCostMaxFlowP(const int src,
                                                    const int snk,
                                                    const vector<int> &ecaps) {

  assert(link_num == 2 * ecaps.size());
  vector<int> dist(vertex_num, INF);
  vector<int> width(vertex_num, 0);
  vector<int> flow(link_num, 0);

  Fixed_heap Q(vertex_num);

  dist[src] = 0;
  width[src] = INF;
  vector<int> caps(link_num);
  for (size_t i = 0; i < ecaps.size(); i++) {
    int link = out2inLink_map[2 * i];
    caps[link] = ecaps[i];
    link = out2inLink_map[2 * i + 1];
    caps[link] = ecaps[i];
  }

  int totflow = 0, totcost = 0;
  int amt = dijkstra(src, snk, caps, dist, width, flow, Q);
  while (amt > 0) {
    totflow += amt;
    for (int x = snk; x != src; x = dad[x].first) {
      int link = (dad[x].second) / 2;

      if (1 == ((dad[x].second) % 2)) {
        flow[link] += amt;
        totcost += amt * link_ends[link].weight;
      } else {
        flow[link] -= amt;
        totcost -= amt * link_ends[link].weight;
      }
    }
    amt = dijkstra(src, snk, caps, dist, width, flow, Q);
  }

  return make_pair(totflow, totcost);
}

void undirected_graph::getMinCostMaxFlow(int src, const int snk,
                                         const vector<int> &ecaps,
                                         vector<vector<int>> &node_paths,
                                         vector<int> &bws,
                                         vector<int> &node_sum_value) {
  assert(link_num == 2 * ecaps.size());
  fill(dist.end(), dist.end(), INF);
  fill(width.begin(), width.end(), 0);
  fill(flow.begin(), flow.end(), 0);
  node_sum_value.resize(vertex_num, 0);
  fill(node_sum_value.begin(), node_sum_value.end(), 0);

  dist[src] = 0;
  width[src] = INF;
  vector<int> caps(link_num);
  for (size_t i = 0; i < ecaps.size(); i++) {
    int link = out2inLink_map[2 * i];
    caps[link] = ecaps[i];
    link = out2inLink_map[2 * i + 1];
    caps[link] = ecaps[i];
  }
  int amt = dijkstra(src, snk, caps);
  while (amt > 0) {
    for (int x = snk; x != src; x = dad[x].first) {
      int link = (dad[x].second) / 2;

      if (1 == ((dad[x].second) % 2)) {
        flow[link] += amt;

      } else {
        flow[link] -= amt;
      }
    }
    amt = dijkstra(src, snk, caps);
  }

  vector<int> path;
  amt = dijkstra(src, snk, flow, path);
  while (amt > 0) {
    int value = 0;
    bws.push_back(amt);
    vector<int> node_path;
    node_path.push_back(src);
    for (vector<int>::iterator it = path.begin(); it != path.end(); it++) {
      int link = *it;
      value += amt * link_ends[link].weight;

      flow[link] -= amt;
      node_path.push_back(link_ends[link].snk);

      node_sum_value[link_ends[link].snk] += value;
    }
    node_paths.push_back(node_path);
    amt = dijkstra(src, snk, flow, path);
  }
}

bool undirected_graph::isValidatePath(const int &src, const int &snk,
                                      const vector<int> &path) const {
  int current = src;
  int next = src;

  for (vector<int>::const_iterator it = path.begin(); it != path.end(); it++) {
    if (!findRhs(2 * (*it), current, next)) {
      return false;
    }

    current = next;
  }

  return next == snk;
}

int undirected_graph::path_cost(const vector<int> &path) const {
  int re = 0;
  for (vector<int>::const_iterator it = path.begin(); it != path.end(); it++) {
    re += link_ends[out2inLink_map[2 * (*it)]].weight;
  }
  return re;
}

void Loc_choose::initial() {

  vector<pair<int, int>> demandC;
  int tempC = totCap;

  for (int v = 0; v < network_node_num; v++) {
    if (network_node_user_map[v] > -1) {
      int user_node = network_node_user_map[v];

      user_to_network_map[user_node] = v;

      int outDegree = graph.getOutDegree(v);
      int sumBw = 0;
      for (int i = 0; i < outDegree; i++) {
        int link = graph.getAdj(v, i);
        sumBw += orignal_caps[link / 2];
      }
      /**
       * if for a directly connect user node when sum of  out link
       *bandwidth less than user demand then this node must be a server
       *
       */

      if (sumBw < 2 * user_demand[user_node]) {
        user_direct_server[user_node] = true;
        choosedServer.insert(v);
        tempC -= user_demand[user_node];
      } else {
        demandC.push_back(make_pair(user_demand[user_node], user_node));
      }
    }
  }

  sort(demandC.rbegin(), demandC.rend());
  for (vector<pair<int, int>>::iterator it = demandC.begin();
       it != demandC.end(); it++) {

    /**
     *  when one user demand greater or equal then sum of all
     * other users demand then this user direct connect node must be
     * a server
     *
     */
    int user_node = it->second;
    int v = user_to_network_map[user_node];

    if (!user_direct_server[user_node]) {
      if (2 * (it->first) >= tempC) {

        user_direct_server[user_node] = true;
        choosedServer.insert(v);
        tempC -= it->first;
      }
    }
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
        break;
      }
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
  for (int v = 0; v < network_node_num; v++) {
    graph.setWeight(2 * v, INF);
    graph.setWeight(2 * v + 1, INF);
  }
  for (int link = target_link_start; link < graph.getLink_num(); link++) {
    graph.setWeight(link, INF);
  }
}
void Loc_choose::releaseVirLinks() {
  for (int v = 0; v < network_node_num; v++) {
    graph.setWeight(2 * v, 0);
    graph.setWeight(2 * v + 1, 0);
  }
  for (int link = target_link_start; link < graph.getLink_num(); link++) {
    graph.setWeight(link, 0);
  }
}
bool Loc_choose::smallestUer() {

  if (2 == user_node_num) {
    int large_user = 0;

    if (user_demand[0] < user_demand[1]) {
      large_user = 1;
    }
    int loc = user_to_network_map[large_user];
    vector<int> caps = orignal_caps;

    caps[loc] = totCap;

    pair<int, int> one_elem =
        graph.getMinCostMaxFlow(virtual_source, virtual_target, caps);

    if (one_elem.first == totCap) {
      Server temp;
      temp.locs.insert(loc);
      temp.success_bw = totCap;
      temp.total_price = one_elem.second + server_price;
      server_candiate.push_back(temp);
    }

    return true;
  }
  if (3 == user_node_num) {
    int large_user = 0;

    if (user_demand[0] < user_demand[1]) {
      large_user = 1;
    }
    if (user_demand[large_user] < user_demand[2]) {
      large_user = 2;
    }

    if (2 * user_demand[large_user] >= totCap) {
    }
  }

  return false;
}
void Loc_choose::tryKServer(const int k) {}

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
  vector<int> node_sum_value;
  graph.getMinCostMaxFlow(virtual_source, virtual_target, caps, node_paths, bws,
                          node_sum_value);

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
char *Loc_choose::solve() {

  /**
   *  there is no user node
   *
   */

  if (0 == user_node_num) {

    char *topo_file = new char[1024];
    fill(topo_file, topo_file + 1023, 0);
    sprintf(topo_file, "0\n\n");
    return topo_file;
  }

  server_candiate.clear();
  Server temp;

  forbidVirLinks();

  for (int user_node = 0; user_node < user_node_num; user_node++) {

    int v = user_to_network_map[user_node];
    temp.locs.insert(v);

    if (user_direct_server[user_node]) {
      continue;
    }

    int limit = server_price / user_demand[user_node];

    vector<pair<int, int>> tree;
    graph.dijkstra_limit_tree(v, limit, tree);
    server_candiate_locs[v].clear();

    for (vector<pair<int, int>>::iterator it = tree.begin(); it != tree.end();
         it++) {
      if (virtual_source != it->first) {
        server_candiate_locs[v].push_back(it->first);
      }
    }

    if (1 == server_candiate_locs[v].size()) {
      user_direct_server[v] = true;
      choosedServer.insert(v);
    }
    sort(server_candiate_locs[v].begin(), server_candiate_locs[v].end());
  }

  releaseVirLinks();

  temp.success_bw = totCap;
  temp.total_price = user_node_num * server_price;
  server_candiate.push_back(temp);

  bool state = true;
  vector<int> haveCheck;

  for (int v = 0; v < user_node_num; v++) {
    std::vector<int> union_data;
    set_union(server_candiate_locs[v].begin(), server_candiate_locs[v].end(),
              haveCheck.begin(), haveCheck.end(),
              std::back_inserter(union_data));
    if (union_data.size() !=
        (server_candiate_locs[v].size() + haveCheck.size())) {
      state = false;
      break;
    }
    sort(union_data.begin(), union_data.end());
    haveCheck = union_data;
  }

  if (state) {
    return output();
  }

  if (smallestUer()) {

    return output();
  }

  value_supper = user_node_num * server_price;

  for (int v = 0; v < network_node_num; v++) {
    vector<int> caps = orignal_caps;
    caps[v] = totCap;
    pair<int, int> one_elem =
        graph.getMinCostMaxFlow(virtual_source, virtual_target, caps);
    Server temp;
    temp.locs.insert(v);
    temp.success_bw = one_elem.first;
    temp.total_price = one_elem.second + server_price;

    server_candiate.push_back(temp);
  }

  for (int v = 0; v < network_node_num; v++) {

    if ((totCap == server_candiate[v].success_bw) &&
        (server_candiate[v].total_price < value_supper)) {
      value_supper = server_candiate[v].total_price;
    }
  }

  /**
   *  one server can support
   *
   */
  if (value_supper <= 2 * server_price) {
    return output();
  }

  value_lower = 2 * server_price;

  for (int k = 0; k < 100; k++) {
    lower_update();
    supper_update();
    if ((value_supper / server_price) == (value_lower / server_price)) {
      return output();
    }
  }

  return output();
}
}

// You need to complete the function
void deploy_server(char *topo[MAX_EDGE_NUM], int line_num, char *filename) {
  int max_line_length = 1024;
  size_t str_len = 1023;
  char *delim = " ";

  char *temp_str = strtok(topo[0], delim);
  int network_node_num = atoi(temp_str);
  temp_str = strtok(NULL, delim);
  int network_link_num = atoi(temp_str);
  temp_str = strtok(NULL, delim);
  int user_node_num = atoi(temp_str);

  int vir_source = network_node_num;
  int vir_target = vir_source + 1;

  vector<int> network_node_user_map(network_node_num, -100);

  temp_str = strtok(topo[2], delim);

  int server_price = atoi(temp_str);

  vector<int> srcs, snks, caps, weights;

  for (int v = 0; v < network_node_num; v++) {
    srcs.push_back(vir_source);
    snks.push_back(v);
    caps.push_back(0);
    weights.push_back(0);
  }

  int line_id = 4;

  int i = 0;
  while (i < network_link_num) {

    i++;
    int src, snk, cap, w;

    temp_str = strtok(topo[line_id++], delim);
    src = atoi(temp_str);
    temp_str = strtok(NULL, delim);
    snk = atoi(temp_str);
    temp_str = strtok(NULL, delim);
    cap = atoi(temp_str);
    temp_str = strtok(NULL, delim);
    w = atoi(temp_str);

    srcs.push_back(src);
    snks.push_back(snk);
    caps.push_back(cap);
    weights.push_back(w);
  }

  line_id++;

  i = 0;
  int target_link_start = 2 * srcs.size();

  vector<int> user_demand;
  while (i < user_node_num) {

    i++;
    int user_node, network_node, bw;
    temp_str = strtok(topo[line_id++], delim);
    user_node = atoi(temp_str);

    temp_str = strtok(NULL, delim);
    network_node = atoi(temp_str);
    temp_str = strtok(NULL, delim);
    bw = atoi(temp_str);

    network_node_user_map[network_node] = user_node;

    srcs.push_back(network_node);
    snks.push_back(vir_target);
    weights.push_back(0);
    caps.push_back(bw);
    user_demand.push_back(bw);
  }

  raptor::undirected_graph graph;
  graph.initial(srcs, snks, weights);

  raptor::Loc_choose loc_choose(
      graph, network_node_num, user_node_num, server_price, vir_source,
      vir_target, target_link_start, network_node_user_map, user_demand, caps);

  char *topo_file = loc_choose.solve();

  write_result(topo_file, filename);
  delete[] topo_file;
}
