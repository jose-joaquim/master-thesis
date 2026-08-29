#include <cstdio>

#include "basic.h"

using namespace std;

int main(int argc, char **argv) {
  freopen(argv[3], "r", stdin);
  read_data();

  freopen(argv[4], "r", stdin);
  int l, ch, ts;
  double tmp;

  vector<pair<int, pair<double, double>>> schedule[512][45];
  while (cin >> l >> ch >> ts >> tmp) {
    schedule[ts][ch].push_back(make_pair(l, make_pair(GMM[l], 0)));
  }

  FILE *log = fopen("log_du8_milp.txt", "a");
  int more = 0;
  for (int i = 0; i < 512; ++i) {
    for (int j = 0; j < 45; ++j) {
      for (int k = 0; k < schedule[i][j].size(); ++k) {
        // printf("conn %d\n", seila[i][j][k].first);
        Connection other_conn(schedule[i][j][k].first);
        for (int u = k + 1; u < schedule[i][j].size(); ++u)
          other_conn.interference +=
              AFF[schedule[i][j][u].first][other_conn.id];

        int bw = idxToB(cToBIdx(j));
        double ans = computeConnectionThroughput(other_conn, bw);
        schedule[i][j][k].second.second = ans;

        // if (!definitelyGreaterThan(ans, GMM[tmp_conn.id]))
        // printf("----- %d %d %d %.3lf %.3lf %.3lf\n", i, j, tmp_conn.id,
        //    tmp_conn.SINR, GMM[tmp_conn.id], ans);

        int minimum_req = 0;
        for (; minimum_req < 4; ++minimum_req)
          if (definitelyLessThan(gammaToBeta(GMM[other_conn.id], minimum_req),
                                 SINR[11][3] + 1.0))
            break;

        assert(minimum_req < 4);
        // if (idxToB(minimum_req) <= bw) {
        //   printf("%d %d %d\n", tmp_conn.id, idxToB(minimum_req), bw);
        // }
        assert(bw >= idxToB(minimum_req));

        fprintf(log, "%d %d %d %d %.3lf %.3lf %s\n", i, other_conn.id,
                idxToB(minimum_req), bw, GMM[other_conn.id], ans, argv[2]);

        // if (idxToB(minimum_req) != bw)
        //   printf("%d %.3lf %d %d %d\n", tmp_conn.id, GMM[tmp_conn.id],
        //          idxToB(minimum_req), j, bw);

        more += idxToB(minimum_req) != bw;
        // definitelyGreaterThan(ans, GMM[tmp_conn.id]);

        assert(idxToB(minimum_req) <= bw);
      }
    }
  }

  fclose(log);
  cout << more << endl;
  return 0;
}
