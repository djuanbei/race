#include "heap.h"

namespace raptor {

void Fixed_heap::percolateUp(int i) {

  
  int w = heap[2*i];
  int loc=heap[2*i+1];
  
  int p = parent(i);
  
  while (i != 0 && (w < heap[2*p])) {
    assgin( i,p );

    indices[heap[2* i+1]] = i;
    i = p;
    p = parent(p);

  }
  
  heap[2*i] = w;
  heap[2*i+1] = loc;
  indices[loc] = i;
  
}

void Fixed_heap::percolateDown(int i) {

  int w=heap[ 2*i ];
  int loc=heap[ 2*i+1 ];

  while (left(i) < size) {
    int child = right(i) < size && (heap[2*right(i)] < heap[2*left(i)])
                    ? right(i)
                    : left(i);
    if (!(heap[2*child] < w)){
      break; 
    }

    assgin( i, child );
    indices[heap[2*i+1]] = i;
    i = child;
  }
  heap[2*i] = w;
  heap[2*i+1] = loc;
  indices[loc] = i;
}

void Fixed_heap::push( int weight, int loc) {
  
  if (-1 == indices[loc]) {
    heap[2*size] = weight;
    heap[2*size+1] = loc;
    
    indices[loc] = size;
    percolateUp(size);
    size++;
  } else {
    
    if (heap[2*indices[loc]] > weight) {
      heap[2*indices[loc]] = weight;
      percolateUp(indices[loc]);
      
    } else if (heap[2*indices[loc]] < weight) {
      heap[2*indices[loc]] = weight;
      percolateDown(indices[loc]);
    }
  }
}

void Fixed_heap::pop() {
  indices[heap[1]] = -1;
  assgin( 0, size-1 );

  if (size > 0) {
    indices[heap[1]] = 0;
    size--;
    percolateDown(0);
  }
}
}
