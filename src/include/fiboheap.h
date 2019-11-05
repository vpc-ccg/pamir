#include <iostream>

using namespace std;

template <class V> class FiboHeap;

template <class V> class Node {
private:
	Node<V>* prev;
	Node<V>* next;
	Node<V>* child;
	Node<V>* parent;
	V value;
	int degree;
	bool marked;
public:
	friend class FiboHeap<V>;
	Node<V>* getPrev() {return prev;}
	Node<V>* getNext() {return next;}
	Node<V>* getChild() {return child;}
	Node<V>* getParent() {return parent;}
	V getValue() {return value;}
	bool isMarked() {return marked;}
	bool hasChildren() {return child;}
	bool hasParent() {return parent;}
};

template <class V> class FiboHeap {
protected:
	// circular doubly linked list
	Node<V>* heap;
	static const int MAXDEGREE = 5000;
public:
	FiboHeap() {
		heap = _empty();
	}
	~FiboHeap() {
		if(heap)
			_deleteAll(heap);
	}
	Node<V>* insert(V value) {
		Node<V>* ret = _singleton(value);
		heap = _merge(heap,ret);
		return ret;
	}
	void merge(FiboHeap& other) {
		heap = _merge(heap,other.heap);
		other.heap = _empty();
	}

	bool isEmpty() {
		return heap == NULL;
	}

	V getMinimum() {
		return heap->value;
	}

	V removeMinimum() {
		Node<V>* old = heap;
		heap = _removeMinimum(heap);
		V ret = old->value;
		delete old;
		return ret;
	}

	void decreaseKey(Node<V>* n,V value) {
		heap = _decreaseKey(heap,n,value);
	}

	void update(Node<V>* n) {
		if (n == heap)
			heap = _getMinIn(heap);
		heap = _update(heap,n);
	}

	V getChild() {
		return heap->child->getValue();
	} 

private:
	Node<V>* _empty() {
		heap = NULL;
		return heap;
	}

	Node<V>* _singleton(V value) {
		Node<V>* n = new Node<V>;
		n->value = value;
		n->prev = n;
		n->next = n;
		n->degree = 0;
		n->marked = false;
		n->child = NULL;
		n->parent = NULL;
		return n;
	}

	// a and b are circular doubly linked lists
	Node<V>* _merge(Node<V>* a,Node<V>* b) {
		if(a == NULL)	return b;
		if(b == NULL)	return a;
		Node<V>* an = a->next;
		Node<V>* bp = b->prev;
		a->next = b;
		b->prev = a;
		an->prev = bp;
		bp->next = an;
		return ((*(a->value) < *(b->value)) ? a : b);
	}

	void _deleteAll(Node<V>* n) {
		if(n == NULL)	return;
		Node<V>* p = n;
		do {
			Node<V>* d = p;
			p = p->next;
			_deleteAll(d->child);
			delete d;
		} while(p != NULL && p != n);
	}

	void _addChild(Node<V>* parent,Node<V>* child) {
		child->prev = child;
		child->next = child;
		child->parent = parent;
		parent->degree++;
		parent->child = _merge(parent->child,child);
	}

	void _unMarkAndUnParentAll(Node<V>* n) {
		if(n == NULL)	return;
		Node<V>* p = n;
		do {
			p->marked = false;
			p->parent = NULL;
			p = p->next;
		}while(p != n);
	}

	Node<V>* _getMinIn(Node<V>* n) {
		Node<V>* min = n;
		Node<V>* p = n;
		do {
			if(*(p->value) < *(min->value))
				min = p;
			p = p->next;
		} while(p != n);
		return min;
	}

	Node<V>* _removeMinimum(Node<V>* n) {
		// meld
		_unMarkAndUnParentAll(n->child);
		if(n->next == n) {
			n = n->child;
		} else {
			n->next->prev = n->prev;
			n->prev->next = n->next;
			n = _merge(n->next,n->child);
		}
		if(n == NULL)	return n;
		
		// consolidate
		Node<V>* trees[MAXDEGREE] = {NULL};
		while(true) {
			if(trees[n->degree] != NULL) {
				Node<V>* t = trees[n->degree];
				if(t == n)	break;
				trees[n->degree] = NULL;
				if(*(n->value) < *(t->value)) {
					t->prev->next = t->next;
					t->next->prev = t->prev;
					_addChild(n,t);
				} else {
					t->prev->next = t->next;
					t->next->prev = t->prev;
					if(n->next == n) {
						t->next = t->prev = t;
						_addChild(t,n);
						n = t;
					} else {
						n->prev->next = t;
						n->next->prev = t;
						t->next = n->next;
						t->prev = n->prev;
						_addChild(t,n);
						n = t;
					}
				}
				continue;
			} else {
				trees[n->degree] = n;
			}
			n = n->next;
		}
		// update min
		return _getMinIn(n);
	}

	Node<V>* _cut(Node<V>* heap,Node<V>* n) {
		if(n->next == n) {
			n->parent->child = NULL;
		} else {
			n->next->prev = n->prev;
			n->prev->next = n->next;
			if (n->parent->child == n)
				n->parent->child = _getMinIn(n->next);
		}
		n->next = n;
		n->prev = n;
		n->marked = false;
		n->parent->degree--;
		n->parent = NULL;
		return _merge(heap,n);
	}

	Node<V>* _decreaseKey(Node<V>* heap, Node<V>* n, V value) {
		if(*(n->value) < *value){
			return heap;
		}
		n->value = value;
		// update structure
		if (n->parent == NULL && *(n->value) < *(heap->value))
			heap = n;
		if(n->parent != NULL && *(n->value) < *((n->parent)->value)) {
			Node<V>* parent = n->parent;
			heap = _cut(heap,n);
			while(parent->parent != NULL && parent->marked) {
				Node<V>* tmp = parent->parent;
				heap = _cut(heap, parent);
				parent = tmp;
			}
			if(parent != NULL && parent->parent != NULL)
				parent->marked = true;
		}
		return heap;
	}
	
	// only appliable for increase key
	Node<V>* _update(Node<V>* heap,Node<V>* n) {
		if (n->child == NULL)	return heap;
		Node<V>* kids = n->child;
		kids = _getMinIn(kids);
		if (*(n->value) < *(kids->value))
			return heap;
		n->child = NULL;
		return _merge(heap, kids);
	}

};
