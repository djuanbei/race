
#ifndef __HEAP_H
#define __HEAP_H

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

}
#endif

