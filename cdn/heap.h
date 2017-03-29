
#ifndef __HEAP_H
#define __HEAP_H

#include <vector>

#include <utility>

namespace raptor {

using namespace std;
class Fixed_heap {
 private:
  vector<int> heap;  // Fixed_heap of Keys
  vector<int> indices;          // Each Key's position (index) in the Fixed_heap

  int size;

  // Index "traversal" functions
  static inline int left(int i) { return (i << 1) | 1; }
  static inline int right(int i) { return (i + 1) << 1; }
  static inline int parent(int i) { return (i - 1) >> 1; }
  inline void assgin(  int i, int j ){
    heap[2*i] = heap[2*j];
    heap[2*i+1] = heap[2*j+1];
  }
  void percolateUp(int i);

  void percolateDown(int i);

 public:
  Fixed_heap(int cap) : heap(2*cap), indices(cap, -1), size(0) {}

  bool empty() const { return 0 == size; }

  void push(int weight, int loc);
  

  void top(int & weight, int & loc) {
    weight=heap[0];
    loc=heap[1];
    // return heap[0];
  }

  // const pair<int, int> &operator[](const int i) const { return heap[i]; }

  void pop();

  int len() const { return size; }

  void clear() {
    size = 0;
    fill(indices.begin(), indices.end(), -1);
  }
};
}
#endif
