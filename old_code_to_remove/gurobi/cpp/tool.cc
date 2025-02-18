#include <cstdio>

#include "header.h"

using namespace std;

int main(int argc, char** argv) {
  freopen(argv[3], "r", stdin);
  read_data();

  // for (int i = 0; i < 12; ++i) {
  //   for (int j = 0; j < 4; ++j) {
  //     printf("%lf ", SINR[i][j]);
  //   }
  //   puts("");
  // }
  //
  // puts("\n");
  //
  // for (int i = 0; i < 12; ++i) {
  //   for (int j = 0; j < 4; ++j) {
  //     printf("%lf ", DR[i][j]);
  //   }
  //   puts("");
  // }

  freopen(argv[4], "r", stdin);
  int l, ch, ts;
  double tmp;

  vector<pair<int, pair<double, double>>> seila[512][45];
  while (cin >> l >> ch >> ts >> tmp) {
    seila[ts][ch].push_back(make_pair(l, make_pair(GMM[l], 0)));
  }

  FILE* log = fopen("log_du8_milp.txt", "a");
  int more = 0;
  for (int i = 0; i < 512; ++i) {
    for (int j = 0; j < 45; ++j) {
      for (int k = 0; k < seila[i][j].size(); ++k) {
        // printf("conn %d\n", seila[i][j][k].first);
        Connection tmp_conn(seila[i][j][k].first);
        for (int u = k + 1; u < seila[i][j].size(); ++u)
          tmp_conn.interference += AFF[seila[i][j][u].first][tmp_conn.id];

        int bw = idxToB(cToBIdx(j));
        double ans = computeConnectionThroughput(tmp_conn, bw);
        seila[i][j][k].second.second = ans;

        // if (!definitelyGreaterThan(ans, GMM[tmp_conn.id]))
        // printf("----- %d %d %d %.3lf %.3lf %.3lf\n", i, j, tmp_conn.id,
        //    tmp_conn.SINR, GMM[tmp_conn.id], ans);

        int minimum_req = 0;
        for (; minimum_req < 4; ++minimum_req)
          if (definitelyLessThan(gammaToBeta(GMM[tmp_conn.id], minimum_req),
                                 SINR[11][3] + 1.0))
            break;

        assert(minimum_req < 4);
        // if (idxToB(minimum_req) <= bw) {
        //   printf("%d %d %d\n", tmp_conn.id, idxToB(minimum_req), bw);
        // }
        assert(bw >= idxToB(minimum_req));

        fprintf(log, "%d %d %d %d %.3lf %.3lf %s\n", i, tmp_conn.id,
                idxToB(minimum_req), bw, GMM[tmp_conn.id], ans, argv[2]);

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
