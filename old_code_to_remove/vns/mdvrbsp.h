#include "common.h"

vector<double> GMA;

bool mdvrbsp_feasible(const Solution &);
optional<Channel> try_insert(const int conn_id, Channel ch);
bool can_split(const Channel &ch);
double computeViolation(const Solution &);

Solution DTS(const Solution &);
Solution RH(Solution);
Solution VNS(Solution);
Solution CH();
