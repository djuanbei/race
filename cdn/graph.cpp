
#include <assert.h>

#include "graph.h"
#include "heap.h"

#include <algorithm>

namespace raptor {
using namespace std;

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
  caps.resize( link_num, 0 );
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

void undirected_graph::dijkstra_tree(const int src, vector<int> &dis) {
  size_t j, outDegree, link, next;
  int current, current_weight;
  int weight;
  dis.resize(vertex_num);
  fill(dis.begin(), dis.end(), INF);
  vector<int> preLink(vertex_num, link_num + 1);
  vector<int> parent(vertex_num, -1);
  dis[src] = 0;
  Fixed_heap Q(vertex_num);

  Q.push(0, src);

  while (!Q.empty()) {
    Q.top( current_weight, current);


    Q.pop();

    outDegree = getOutDegree(current);
    // int current_weight = p.first;
    for (j = 0; j < outDegree; j++) {
      link = outIndex[current] + j;
      const endElement &neighbour = link_ends[link];
      weight = current_weight + neighbour.weight;
      next = neighbour.snk;
      if (weight < dis[next]) {
        parent[next] = current;
        preLink[next] = link;
        dis[next] = weight;
        Q.push(weight, next);
      }
    }
  }
}

void undirected_graph::dijkstra_limit_tree(const int src, const int limit,
                                           vector<pair<int, int>> &tree) {
  tree.clear();
  if (0 == limit) {
    return;
  }

  size_t j, outDegree, link, next;
  int current;
  int weight, current_weight;
  vector<int> dis(vertex_num, INF);
  vector<int> preLink(vertex_num, link_num + 1);
  vector<int> parent(vertex_num, -1);
  dis[src] = 0;
  Fixed_heap Q(vertex_num);

  Q.push(0, src);

  while (!Q.empty()) {
    Q.top( current_weight, current);

    Q.pop();

    outDegree = getOutDegree(current);

    for (j = 0; j < outDegree; j++) {
      link = outIndex[current] + j;
      const endElement &neighbour = link_ends[link];
      weight = current_weight + neighbour.weight;
      next = neighbour.snk;
      if (weight <= limit && weight < dis[next]) {
        parent[next] = current;
        preLink[next] = link;
        dis[next] = weight;
        Q.push(weight, next);
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
    Q.push(-outDegree, i);
    degrees[i] = -outDegree;
  }

  vector<bool> check(link_num, false);

  while (!Q.empty()) {
    Q.top(next, current);
    
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
        Q.push(degrees[next], next);
      }
    }
  }
}

int undirected_graph::dijkstra(const int src, const int snk,
                               const vector<int> &caps, vector<int> &dist,
                               vector<int> &width, vector<int> &flow,
                               Fixed_heap &Q) {
  int j, outDegree, link, next;
  int current;
  int weight, current_weight;

  fill(dist.begin(), dist.end(), INF);
  fill(width.begin(), width.end(), 0);
  dist[src] = 0;
  width[src] = INF;

  Q.clear();
  Q.push(0, src);

  while (!Q.empty()) {
    Q.top( current_weight, current);

    if (current == snk) {
      return width[snk];
    }
    Q.pop();

    outDegree = getOutDegree(current);



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
          Q.push(weight, next);
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
          Q.push(weight, next);
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
  int weight, current_weight;
  vector<int> dis(vertex_num, INF);
  vector<int> preLink(vertex_num, link_num + 1);
  vector<int> parent(vertex_num, -1);
  dis[src] = 0;
  Q.clear(  );
  
  // Fixed_heap Q(vertex_num);

  Q.push(0, src);

  while (!Q.empty()) {
    Q.top( current_weight, current);

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
          Q.push(weight, next);
        }
      }
    }
  }

  return 0;
}

pair<int, int> undirected_graph::getMinCostMaxFlow(
    const int src, const int snk, const vector<int> &ecaps,  const int leftCap,
    vector<int> &outputs, vector<int> &inputs, vector<int> &node_sum_value, vector<vector<int> > &sum_of_pass_flow) {
  
  
  fill(flow.begin(), flow.end(), 0);
  outputs.resize(vertex_num, 0);
  fill(outputs.begin(), outputs.end(), 0);

  inputs.resize(vertex_num, 0);
  fill(inputs.begin(), inputs.end(), 0);

  node_sum_value.resize(vertex_num, 0);
  fill(node_sum_value.begin(), node_sum_value.end(), 0);
  
  // vector<int> caps(link_num);
  for (size_t i = 0; i < ecaps.size(); i++) {
    int link = out2inLink_map[2 * i];
    caps[link] = ecaps[i];
    link = out2inLink_map[2 * i + 1];
    caps[link] = ecaps[i];
  }

  int totflow = 0, totcost = 0;
  int amt = dijkstra(src, snk, caps);
  while (amt > 0) {
    
    int lastLink = (dad[snk].second) / 2;
    inputs[_srcs[lastLink]] += amt;

    int link=0;
    
    totflow += amt;
    
    for (int x = snk; x != src; x = dad[x].first) {
      
      link = (dad[x].second) / 2;

      if (1 == ((dad[x].second) % 2)) {
        flow[link] += amt;
        totcost += amt * link_ends[link].weight;
      } else {
        flow[link] -= amt;
        totcost -= amt * link_ends[link].weight;
      }
    }

    outputs[link_ends[link].snk] += amt;
    
    amt = dijkstra(src, snk, caps);
  }
  
  if(leftCap== totflow ){
    
    vector<int> path;
    amt = dijkstra(src, snk, flow, path);
    while (amt > 0) {
      int network_node_server=link_ends[path.front(  )].snk;
      int value = 0;

      for (vector<int>::iterator it = path.begin()+1; it+1 != path.end(); it++) {
        int link = *it;
        value += amt * link_ends[link].weight;

        flow[link] -= amt;

        node_sum_value[link_ends[link].snk] += value;

        sum_of_pass_flow[ network_node_server ][link_ends[link].snk  ]+=amt;
      }

      amt = dijkstra(src, snk, flow, path);
    }
  }
  
  return make_pair(totflow, totcost);
}

void undirected_graph::getMinCostMaxFlow(int src, const int snk,
                                         const vector<int> &ecaps,
                                         vector<vector<int>> &node_paths,
                                         vector<int> &bws) {
  assert(link_num == 2 * ecaps.size());
  fill(dist.end(), dist.end(), INF);
  fill(width.begin(), width.end(), 0);
  fill(flow.begin(), flow.end(), 0);

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
    bws.push_back(amt);
    vector<int> node_path;
    node_path.push_back(src);
    for (vector<int>::iterator it = path.begin(); it != path.end(); it++) {
      int link = *it;

      flow[link] -= amt;
      node_path.push_back(link_ends[link].snk);
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
}
