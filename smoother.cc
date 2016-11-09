/// 786

#include <bits/stdc++.h>
#include <fmt/format.h>
#include <boost/heap/fibonacci_heap.hpp>
using namespace std;
using namespace boost::heap;

#define E(msg,...)\
	fmt::print(stderr, msg"\n", ##__VA_ARGS__)
#define L(msg,...)\
	fmt::print(stdout, msg"\n", ##__VA_ARGS__)

unordered_map<string, int> read_to_id;
struct Xread {
	vector<int> clusters;
};
struct cluster {
	int id;
	int orig_id;
	string name;
	int support, unresolved;
	vector<int> reads;
};

struct cluster_compare {
	bool operator() (const cluster *a, const cluster *b) const {
		return make_pair(a->support, a->unresolved) < make_pair(b->support, b->unresolved);
	}
};

void set_cover (auto &clusters, auto &reads) 
{
	fibonacci_heap<cluster*, compare<cluster_compare>> heap;
	unordered_map<int, decltype(heap)::handle_type> heap_handles;

	unordered_map<int,bool> x;
	for (auto &c: clusters) {
		if (c.unresolved == 0) continue;
		heap_handles[c.id] = heap.push(&c);
		x[c.id] = false;
	}

	// assert that all shared clusters are in set cover
	int re = 0;
	for (auto &r: reads) {
		assert(r.clusters.size() >= 1);
		if (r.clusters.size() > 1) {
			re++;
			for (auto &c: r.clusters) {
				assert(heap_handles.find(c) != heap_handles.end());
			}
		}
	}

	unordered_set<int> W;

	while (!heap.empty()) {
		cluster *c = heap.top(); heap.pop();
		if (!c->unresolved && !c->support) break;
		//E("{}-- {} {}", c->id, c->support, c->unresolved);

		c->unresolved = 0; // we will resolve all reads!
		heap_handles.erase(c->id);
		x[c->id] = true;

		//if (!c->unresolved) continue;

		for (auto &r: c->reads) {
			if (reads[r].clusters.size() <= 1) continue;

			for (auto &cr: reads[r].clusters) {
				auto it = heap_handles.find(cr);
				if (it == heap_handles.end()) continue; // not in heap anymore
				assert((*it->second)->id == cr);
				assert((*it->second)->unresolved > 0);

				(*it->second)->unresolved--;
				heap.update(it->second);
			}
			reads[r].clusters.clear();
			reads[r].clusters.push_back(c->id);
		}
	}

	for (auto &c: clusters) {
		c.reads.clear();
	}

	unordered_map<int, string> read_names;
	for (auto &r: read_to_id)
		read_names[r.second] = r.first;
	for (auto &r: reads) {
		assert(r.clusters.size() == 1);
		clusters[r.clusters[0]].reads.push_back(&r - &reads[0]);
	}
	int i = 0, rx = 0;
	for (auto &c: clusters) {
		if (c.reads.size() == 0) {
			L("Removed: {} {}", c.orig_id, c.name);
			continue;
		}

		L("{} {} {}", c.orig_id, c.reads.size(), c.name);
		for (auto &r: c.reads) {
			//L("{}", read_names[r]);
		}
		
		rx += c.reads.size(), i++;
	}
	assert(rx == reads.size());
	E("{} / {} sets resolved, {} / {} reads resolved", i, clusters.size(), re, rx);
}

int main (int argc, char **argv) 
{
	freopen(argv[1], "r", stdin);

	vector<cluster> clusters;
	vector<Xread> reads;

	FILE *fi = fopen(argv[1], "r");
	fseek(fi, 0, SEEK_END);
	double fsz = ftell(fi);
	fseek(fi, 0, 0);

	int cid, num, st, ed;
	char name[500];
	char clid[500];
	while (fscanf(fi, "%d %d %s", &cid, &num, clid) != EOF) {
		int id = clusters.size();
		clusters.push_back({id, cid, string(clid), 0, 0, vector<int>()});
		for (int i = 0; i < num; i++) {
			fscanf(fi, "%s", name);
			int rid = read_to_id.size();
			auto f = read_to_id.find(string(name));
			if (f == read_to_id.end()) {
				reads.push_back({vector<int>()});
				read_to_id[name] = rid;
			} else {
				rid = f->second;
			}
			reads[rid].clusters.push_back(id);
			clusters[id].reads.push_back(rid);
		}

		if (clusters.size() % 1000 == 1) {
			auto p = 100 * ftell(fi) / fsz;
			fmt::print(stderr, "\r {:.2f}", p);
		}
	}
	E("");

	for (auto &c: clusters) {
		auto s = unordered_set<int>(c.reads.begin(), c.reads.end());
		c.reads = vector<int>(s.begin(), s.end());
		for (auto &r: c.reads) {
			auto s = unordered_set<int>(reads[r].clusters.begin(), reads[r].clusters.end());
			reads[r].clusters = vector<int>(s.begin(), s.end());
			assert(reads[r].clusters.size() > 0);
			c.unresolved += reads[r].clusters.size() > 1;
		}
		c.support = c.reads.size() - c.unresolved;
	}

	E("{} clusters", clusters.size());
	E("{} reads", reads.size());

	set_cover(clusters, reads);

	// for (auto &r: reads) {
	// 	if (r.clusters.size() != 1) {
	// 		E("{} {}", (&r - &reads[0]), r.clusters.size());
	// 	}
	// }

	E("done");

	return 0;
}
