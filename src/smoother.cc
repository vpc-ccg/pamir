/// 786

#include "fiboheap.h"
#include <cassert>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "logger.h"

using namespace std;
namespace smoother {
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

        bool operator<(const cluster &other) const {
            //return make_pair(support, unresolved) > make_pair(other.support, other.unresolved);
            return make_pair(support, make_pair(unresolved, id)) >
                   make_pair(other.support, make_pair(other.unresolved, other.id));
        }

        friend ostream &operator<<(ostream &os, const cluster &cl);
    };

    ostream &operator<<(ostream &os, const cluster &cl) {
        os << cl.id << '\t' << cl.orig_id << '\t' << cl.name << '\t' << cl.support << '\t' << cl.unresolved << '\t'
           << cl.reads.size();
        return os;
    }

//void set_cover (auto &clusters, auto &reads) 
    void set_cover(vector <cluster> &clusters, vector <Xread> &reads) {
        FiboHeap < cluster * > heap;
        unordered_map < int, Node < cluster * > * > heap_handles;

        // inserting them all
        unordered_map<int, bool> x;
        for (auto &c: clusters) {
            if (c.unresolved == 0) continue;
            heap_handles[c.id] = heap.insert(&c);
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

        while (!heap.isEmpty()) {
            cluster *c = heap.getMinimum();
            heap.removeMinimum();
            //cout << c->id << '\t' << c->support <<'\t' << c->unresolved << '\t' << c->reads.size()<<'\n';

            if (!c->unresolved && !c->support) break;
            //E("{}-- {} {}", c->id, c->support, c->unresolved);

            c->unresolved = 0; // we will resolve all reads!
            heap_handles.erase(c->id);
            x[c->id] = true;

            //if (!c->unresolved) continue;

            for (auto &r: c->reads) {
                if (reads[r].clusters.size() <= 1) continue;

                for (auto &cr: reads[r].clusters) {
                    unordered_map < int, Node < cluster * > * > ::const_iterator
                    it = heap_handles.find(cr);
                    if (it == heap_handles.end()) continue; // not in heap anymore

                    assert(((it->second)->getValue())->id == cr);
                    assert(((it->second)->getValue())->unresolved > 0);

                    cluster *newc = it->second->getValue();
                    newc->unresolved--;
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
                Logger::instance().info("Removed: %d %s\n", c.orig_id, c.name.c_str());
                continue;
            }

            Logger::instance().info("%d %lu %s\n", c.orig_id, c.reads.size(), c.name.c_str());
            for (auto &r: c.reads) {
                //L("{}", read_names[r]);
            }

            rx += c.reads.size(), i++;
        }
        assert(rx == reads.size());
        Logger::instance().error("%d / %lu sets resolved, %d / %d reads resolved\n", i, clusters.size(), re, rx);
    }

    int main(int argc, char **argv) {
        freopen(argv[1], "r", stdin);

        vector <cluster> clusters;
        vector <Xread> reads;

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
                Logger::instance().error("\r %f\n", p);
            }
        }
        Logger::instance().error("\n");

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

        Logger::instance().error("%lu clusters\n", clusters.size());
        Logger::instance().error("%lu reads\n", reads.size());

        set_cover(clusters, reads);

        // for (auto &r: reads) {
        // 	if (r.clusters.size() != 1) {
        // 		E("{} {}", (&r - &reads[0]), r.clusters.size());
        // 	}
        // }

        Logger::instance().error("done\n");

        return 0;
    }
}