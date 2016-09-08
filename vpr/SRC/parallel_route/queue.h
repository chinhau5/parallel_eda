#ifndef QUEUE_H
#define QUEUE_H

template<typename T>
struct queue_t {
	T *items;
	int capacity;
	int read;
	int write;
};

template<typename T>
void q_init(queue_t<T> *q, int capacity)
{
	q->capacity = capacity+1;
	q->items = new T[q->capacity];
	q->read = 0;
	q->write = 0;
}

template<typename T>
int q_size(const queue_t<T> *q)
{
	return (q->capacity + q->write - q->read) % q->capacity;
}

template<typename T>
bool q_full(const queue_t<T> *q)
{
	return q->read == (q->write+1) % q->capacity;
}

template<typename T>
bool q_empty(const queue_t<T> *q)
{
	return q->read == q->write;
}

template<typename T>
T &q_front(queue_t<T> *q)
{
	assert(!q_empty(q));
	return q->items[q->read];
}

template<typename T>
T &q_back(queue_t<T> *q)
{
	assert(!q_empty(q));
	return q->items[(q->capacity + q->write - 1) % q->capacity];
}

template<typename T>
bool q_push(queue_t<T> *q, const T &item)
{
	if (q_full(q)) {
		return false;
	}
	q->items[q->write] = item;
	q->write = (q->write+1) % q->capacity;
	return true;
}

template<typename T>
bool q_pop(queue_t<T> *q)
{
	if (q_empty(q)) {
		return false;
	}
	q->read = (q->read+1) % q->capacity;
	return true;
}

#endif
