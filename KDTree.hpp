// Copyright

#ifndef SRC_KDTREE_HPP_
#define SRC_KDTREE_HPP_

#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>
#include "Point.hpp"

using namespace std;

template <size_t N, typename ElemType> 
struct KDTreeNode {
  KDTreeNode<N,ElemType>* a_nodes[2];
  ElemType value;
  Point<N> pt_here;
  KDTreeNode(const Point<N> p) {
    a_nodes[0] = nullptr;
    a_nodes[1] = nullptr;
    pt_here = p;
  }
  KDTreeNode(const Point<N> p, const ElemType el) {
    a_nodes[0] = nullptr;
    a_nodes[1] = nullptr;
    pt_here = p;
    value = el;
  }
};

template <size_t N, typename ElemType>
class KDTree {
 public:
  typedef std::pair<Point<N>, ElemType> value_type;

  KDTree();

  ~KDTree();

  KDTree(const KDTree &rhs);
  KDTree &operator=(const KDTree &rhs);

  size_t dimension() const;

  size_t size() const;
  bool empty() const;

  bool find(KDTreeNode<N,ElemType>** &p, const Point<N> &pt);

  const bool find2 (KDTreeNode<N,ElemType>*const* &p, const Point<N> &pt) const;

  bool contains(const Point<N> &pt) const;

  void insert(const Point<N> &pt, const ElemType &value);

  ElemType &operator[](const Point<N> &pt);

  ElemType &at(const Point<N> &pt);
  const ElemType &at(const Point<N> &pt) const;

  ElemType knn_value(const Point<N> &key, size_t k) const;

  std::vector<ElemType> knn_query(const Point<N> &key, size_t k) const;

 private:
  KDTreeNode<N,ElemType>* root;
  size_t dimension_;
  size_t size_;
};

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree() {
  root = nullptr;
  dimension_ = N;
  size_ = 0;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::~KDTree() {
  // TODO(me): Fill this in.
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree(const KDTree &rhs) {
  // TODO(me): Fill this in.
}

template <size_t N, typename ElemType>
KDTree<N, ElemType> &KDTree<N, ElemType>::operator=(const KDTree &rhs) {
  return *this;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::dimension() const {
  return dimension_;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::size() const {
  return size_;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::empty() const {
  return !root;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::find(KDTreeNode<N,ElemType>** &p, const Point<N> &pt) {
  int axis = 0;
  for (p = &root; *p && ((*(*p)).pt_here != pt); p = &((*p)->a_nodes[(*(*p)).pt_here[axis] < pt[axis]]), axis++) {//Recorre el arbol hasta que el puntero apunte a nulo, o se encuentre el valor
    //cout << "wtf?\naxis -> " << axis << endl;
    if (axis == N-1) {
      axis = 0;
    }
  }
	return !!*p;
}

template <size_t N, typename ElemType>
const bool KDTree<N, ElemType>::find2(KDTreeNode<N,ElemType>*const* &p, const Point<N> &pt) const {
  int axis = 0;
  for (p = &root; *p && ((*(*p)).pt_here != pt); p = &((*p)->a_nodes[(*(*p)).pt_here[axis] < pt[axis]]), axis++) {
    if (axis == N-1) {
      axis = 0;
    }
  }
	return !!*p;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::contains(const Point<N> &pt) const {
  KDTreeNode<N,ElemType>*const* p = nullptr;
  if (find2(p, pt)) { return true; }
  return false;
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::insert(const Point<N> &pt, const ElemType &value) {
  KDTreeNode<N,ElemType>** p = nullptr;
	if (find(p, pt)) { (*p)->value = value;; return; }
	*p = new KDTreeNode<N,ElemType> (pt, value);
  size_++;
}

template <size_t N, typename ElemType>
ElemType &KDTree<N, ElemType>::operator[](const Point<N> &pt) {
  KDTreeNode<N,ElemType>** p = nullptr;
	if (find(p, pt)) { return (*p)->value; }
	*p = new KDTreeNode<N,ElemType> (pt);
  size_++;
  return (*p)->value;
}

template <size_t N, typename ElemType>
ElemType &KDTree<N, ElemType>::at(const Point<N> &pt) {
  KDTreeNode<N,ElemType>** p = nullptr;
  if (find(p, pt)) { return (*(*p)).value; }
  throw out_of_range("The element does not exist.\n");
}

template <size_t N, typename ElemType>
const ElemType &KDTree<N, ElemType>::at(const Point<N> &pt) const {
  KDTreeNode<N,ElemType>*const* p = nullptr;
  if (find2(p, pt)) { return (*(*p)).value; }
  throw out_of_range("The element does not exist.\n");
}

template <size_t N, typename ElemType>
ElemType KDTree<N, ElemType>::knn_value(const Point<N> &key, size_t k) const {
  // TODO(me): Fill this in.
  ElemType new_element;
  return new_element;
}

template <size_t N, typename ElemType>
std::vector<ElemType> KDTree<N, ElemType>::knn_query(const Point<N> &key,
                                                     size_t k) const {
  // TODO(me): Fill this in.
  std::vector<ElemType> values;
  return values;
}

// TODO(me): finish the implementation of the rest of the KDTree class

#endif  // SRC_KDTREE_HPP_