#include "heap.h"

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
    if (!(heap[child].first < x.first)) break;
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
    
    if(heap[indices[k.second]].first> k.first  ){
      heap[indices[k.second]].first = k.first;
      percolateUp(indices[k.second]);
    }else if( heap[indices[k.second]].first < k.first  ){
      heap[indices[k.second]].first = k.first;
      percolateDown(indices[k.second]);
    }
    
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


}
