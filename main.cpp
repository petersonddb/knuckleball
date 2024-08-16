#include <fstream>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

const double deg2rad = M_PI / 180;
const double rev2deg = 360;

using json = nlohmann::json;

struct vector3 {
  double x, y, z;

  vector3 operator*(double scale);
  vector3 operator+(vector3 op);

  double length();
};

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(vector3, x, y, z)

vector3 vector3::operator*(double scale) {
  return {x * scale, y * scale, z * scale};
}

vector3 vector3::operator+(vector3 op) {
  return {x + op.x, y + op.y, z + op.z};
}

double vector3::length() { return sqrt(x * x + y * y + z * z); }

struct params {
  vector3 v0;
  double wy, theta0, xt, time_step;
  const vector3 x0;
  const double g, s0m;

  params() : g(9.81), s0m(4.1e-4), x0({0, 1.1, 0}){};
};

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(params, v0, wy, theta0, xt, time_step)

struct data {
  std::vector<vector3> x, v;
};

void initialize(params &p, const std::string file_name);
void calculate(data &d, const params &p);
void store_data(const data &d, const std::string file_name);

int main(int argc, char **argv) {
  params p;
  data d;

  initialize(p, std::string(argv[1]));
  calculate(d, p);
  store_data(d, std::string(argv[1]).append(".dat"));
}

void initialize(params &p, std::string file_name) {
  std::ifstream file(file_name);
  json read_json{json::parse(file)};
  read_json.get_to(p);
}

void calculate(data &d, const params &p) {
  double theta = p.theta0;
  d.x.push_back(p.x0);
  d.v.push_back(p.v0);
  while (d.x.back().x < p.xt) {
    d.x.push_back(d.x.back() + d.v.back() * p.time_step);

    double b2m = 0.0039 + 0.0058 / (1 + exp((d.v.back().length() - 35) / 5));

    double flm =
        p.g * 0.5 *
        (sin(4 * theta * deg2rad) - 0.25 * sin(8 * theta * deg2rad) +
         0.08 * sin(12 * theta * deg2rad) - 0.0225 * sin(16 * theta * deg2rad));

    d.v.push_back(d.v.back() +
                  vector3{-b2m * d.v.back().length() * d.v.back().x, -p.g,
                          p.s0m * d.v.back().x * p.wy + flm} *
                      p.time_step);

    theta = theta + p.wy * rev2deg * p.time_step;
  }
}

void store_data(const data &d, const std::string file_name) {
  std::ofstream file(file_name);
  int n = d.x.size();

  for (int i = 0; i < n; i++) {
    file << d.x[i].x << " " << d.x[i].y << " " << d.x[i].z << " " << d.v[i].x
         << " " << d.v[i].y << " " << d.v[i].z << std::endl;
  }
}
