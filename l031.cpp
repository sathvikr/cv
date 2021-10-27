#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <list>
#include <cfloat>
#include <unordered_set>
#include <queue>
#include <fstream>
#include <algorithm>

using namespace std;

class Point {

private:
	double x;
	double y;

public:
	Point() {
		x = 0;
		y = 0;
	}

	Point(double x, double y) {
		this->x = x;
		this->y = y;
	}

	double get_x() {
		return x;
	}

	void set_x(int x) {
		this->x = x;
	}

	double get_y() {
		return y;
	}

	void set_y(int y) {
		this->y = y;
	}

	double get_distance_from(Point other) {
		return sqrt(
				pow(other.get_y() - this->get_y(), 2)
						+ pow(other.get_x() - this->get_x(), 2));
	}

	bool equals(Point other) {
		return other.get_x() == x && other.get_y() == y;
	}

	bool imprecise_equals(Point other) {
		return get_display_string() == other.get_display_string();
	}

	Point copy() {
		return Point(x, y);
	}

	string get_display_string() {
		return "(" + to_string(x) + ", " + to_string(y) + ")";
	}

	void print() {
		cout << get_display_string() << endl;
	}

};

class Line {

private:
	Point p1;
	Point p2;

public:
	Line(Point p1, Point p2) {
		this->p1 = p1;
		this->p2 = p2;
	}

	Point get_p1() {
		return p1;
	}

	void set_p1(Point p1) {
		this->p1 = p1;
	}

	Point get_p2() {
		return p2;
	}

	void set_p2(Point p2) {
		this->p2 = p2;
	}

	double get_slope() {
		if (p2.get_x() - p1.get_x() == 0) {
			return DBL_MAX;
		}

		return ((double) p2.get_y() - p1.get_y())
				/ ((double) p2.get_x() - p1.get_x());
	}

	double get_y_intercept() {
		return -get_slope() * p1.get_x() + p1.get_y();
	}

	double get_perpendicular_slope() {
		double m = get_slope();

		if (m == 0) {
			return DBL_MAX;
		} else if (m == DBL_MAX) {
			return 0;
		} else {
			return -1 / m;
		}
	}

	double get_length() {
		return p1.get_distance_from(p2);
	}

	Point get_point_at(double x) {
		return Point(x, get_slope() * x + get_y_intercept());
	}

	Point get_shifted_point(double distance) {
		double x1 = p1.get_x();
		double x2 = p2.get_x();
		double y1 = p1.get_y();
		double y2 = p2.get_y();

		double x3 = distance * (x2 - x1) / get_length() + x1;
		double y3 = (y2 - y1) * (x3 - x1) / (x2 - x1) + y1;

		return Point(x3, y3);
	}

	Point get_intersection(Line other) {
		double a1 = (double) (this->get_p2().get_y() - this->get_p1().get_y());
		double b1 = (double) (this->get_p1().get_x() - this->get_p2().get_x());
		double c1 = a1 * (this->get_p1().get_x())
				+ b1 * (this->get_p1().get_y());

		double a2 = (double) (other.get_p2().get_y() - other.get_p1().get_y());
		double b2 = (double) (other.get_p1().get_x() - other.get_p2().get_x());
		double c2 = a2 * (other.get_p1().get_x())
				+ b2 * (other.get_p1().get_y());

		double determinant = a1 * b2 - a2 * b1;

		double x = (b2 * c1 - b1 * c2) / determinant;
		double y = (a1 * c2 - a2 * c1) / determinant;

		return Point(x, y);
	}

	Line get_perpendicular(Point p) {
		double perpendicular_slope = get_perpendicular_slope();

		if (perpendicular_slope == DBL_MAX) {
			return Line(Point(p.get_x(), 0), Point(p.get_x(), 1));
		} else {
			double b = -perpendicular_slope * p.get_x() + p.get_y();

			return Line(Point(0, b), Point(1, perpendicular_slope * 1 + b));
		}
	}

	void extend(double x1, double x2) {
		this->set_p1(get_point_at(x1));
		this->set_p2(get_point_at(x2));
	}

	bool equals(Line other) {
		return (p1.equals(other.get_p1()) || p1.equals(other.get_p2()))
				&& (p2.equals(other.get_p1()) || p2.equals(other.get_p2()));
	}

	bool imprecise_equals(Line other) {
		return p1.imprecise_equals(other.get_p1())
				&& p2.imprecise_equals(other.get_p2());
	}

	string get_display_string() {
		return p1.get_display_string() + ", " + p2.get_display_string();
	}

};

int pixels[800][800][3];
int size = (int) sizeof(pixels) / sizeof(pixels[0]);

void set_pixel(int x, int y, int r, int g, int b) {
	if (x >= 0 && y >= 0 && x < 800 && y < 800) {
		pixels[x][y][0] = r;
		pixels[x][y][1] = g;
		pixels[x][y][2] = b;
	}
}

void draw_circle(int a, int b, double r, int red, int green, int blue) {
	double radius = r;
	int x, y, xmax, y2, y2_new, ty;

	xmax = (int) (radius * 0.70710678);
	y = r;
	y2 = y * y;
	ty = (2 * y) - 1;
	y2_new = y2;

	for (x = 0; x <= xmax; x++) {
		if ((y2 - y2_new) >= ty) {
			y2 -= ty;
			y -= 1;
			ty -= 2;
		}

		set_pixel(a + x, b + y, red, green, blue);
		set_pixel(a + x, b - y, red, green, blue);
		set_pixel(a - x, b + y, red, green, blue);
		set_pixel(a - x, b - y, red, green, blue);
		set_pixel(a + y, b + x, red, green, blue);
		set_pixel(a + y, b - x, red, green, blue);
		set_pixel(a - y, b + x, red, green, blue);
		set_pixel(a - y, b - x, red, green, blue);

		y2_new -= (2 * x) - 3;
	}
}

void plot_point(Point p) {
	draw_circle(p.get_x() * size, p.get_y() * size, 3, 0, 0, 0);
}

void plot_highlighted_point(Point p, int r, int g, int b) {
	draw_circle(p.get_x() * size, p.get_y() * size, 2, r, g, b);
	draw_circle(p.get_x() * size, p.get_y() * size, 3, r, g, b);
}

void clear_pixels() {
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			for (int k = 0; k < 3; k++) {
				pixels[i][j][k] = 255;
			}
		}
	}
}

void write_to_ppm(string filename) {
	ofstream outfile(filename, ios_base::out | ios_base::binary);
	outfile << "P3\n";
	outfile << size << " " << size << endl;
	outfile << "255\n";

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			for (int k = 0; k < 3; k++) {
				outfile << pixels[i][j][k] << " ";
			}
		}

		outfile << "\n";
	}

	outfile.close();
}

double random_double() {
	static bool first = true;

	if (first) {
		srand(time(NULL));
		first = false;
	}

	return (double) rand() / RAND_MAX;
}

Point get_random_point() {
	return Point(random_double(), random_double());
}

list<Point> get_n_points_list(int n) {
	list<Point> points;

	int i = 0;

	while (i < n) {
		points.push_back(get_random_point());

		i++;
	}

	return points;
}

vector<Point> get_n_points_vector(int n) {
	vector<Point> points;

	int i = 0;

	while (i < n) {
		points.push_back(get_random_point());

		i++;
	}

	return points;
}

list<Point> get_closest_points_bruteforce(list<Point> points) {
	list<Point> closest_points;
	Point closest1, closest2;
	double min_distance = DBL_MAX;

	unordered_set<string> visited;

	list<Point>::iterator p1;
	list<Point>::iterator p2;

	for (p1 = points.begin(); p1 != points.end(); ++p1) {
		int distance = (int) std::distance<list<Point>::iterator>(points.begin(), p1);

		for (p2 = next(points.begin(), distance + 1); p2 != points.end(); ++p2) {
			Line curr = Line(*p1, *p2);
			Line reversed = Line(*p2, *p1);

			if (visited.size() == 0
					|| visited.count(reversed.get_display_string()) == 0) {
				double distance = curr.get_length();

				if (distance < min_distance) {
					closest1 = *p1;
					closest2 = *p2;

					min_distance = distance;
				}

				visited.insert(reversed.get_display_string());
			}
		}
	}

	closest_points.push_back(closest1);
	closest_points.push_back(closest2);

	cout << min_distance << endl;

	return closest_points;
}

double get_min_dist_bruteforce(vector<Point> points) {
	list<Point> closest_points;
	double min_distance = DBL_MAX;

	unordered_set<string> visited;

	for (int i = 0; i < points.size(); i++) {
		for (int j = i + 1; j < points.size(); j++) {
			Line curr = Line(points[i], points[j]);
			Line reversed = Line(points[j], points[i]);

			if (visited.size() == 0
					|| visited.count(reversed.get_display_string()) == 0) {
				double distance = curr.get_length();

				if (distance < min_distance) {
					min_distance = distance;
				}

				visited.insert(reversed.get_display_string());
			}
		}
	}

	return min_distance;
}

bool comparator(Point& lhs, Point& rhs) {
   return lhs.get_x() < rhs.get_x();
}

void write_and_plot(string filename, list<Point> points, int precision) {
	std::ofstream f;

	f.open(filename);
	f << fixed << setprecision(precision);

	for (Point p : points) {
		plot_point(p);
		f << p.get_x() << "  " << p.get_y() << endl;
	}

	f.close();
}

vector<Point> get_strip(vector<Point> points, double half, double width) {
	vector<Point> strip;
	Point midpoint = points[half];

	for (int i = 0; i < points.size(); i++) {
		if (abs(points[i].get_x() - midpoint.get_x()) < width) {
			strip.push_back(points[i]);
		}
	}

	return strip;
}

double get_min_dist_strip(vector<Point> strip) {
	double min_dist = DBL_MAX;
	double half = strip.size() / 2;

	for (int i = 0; i < half; i++) {
		for (int j = half; j < strip.size(); j++) {
			double curr_dist = strip[i].get_distance_from(strip[j]);

			if (curr_dist < min_dist) {
				min_dist = curr_dist;
			}
		}
	}

	return min_dist;
}


double get_min_dist_recur(int start, int end, vector<Point> points) {

	if (end - start == 2) {
		return points[start].get_distance_from(points[start + 1]);
	} else if (end - start == 3) {
		vector<Point> temp;

		temp.push_back(points[start]);
		temp.push_back(points[start + 1]);
		temp.push_back(points[start + 2]);

		return get_min_dist_bruteforce(temp);
	}

	double half = (start + end) / 2;

	double left_min = get_min_dist_recur(start, half, points); // [0, 6)
	double right_min = get_min_dist_recur(half, end, points); // [6, 12)


//	cout << "mins: " << left_min << ", " << right_min << endl;

	double min_dist = min(left_min, right_min);

//	cout << min_dist << endl;

	vector<Point> strip = get_strip(points, half, min_dist);

	double min_strip_dist = get_min_dist_bruteforce(strip);

	return min(min_dist, min_strip_dist);
}

std::string remove_characters(std::string str, char c) {
	str.erase(std::remove(str.begin(), str.end(), c), str.end());

    return str;
}

std::vector<std::string> split(std::string s, std::string del) {
    std::vector<std::string> splitted;

    int start = 0;
    int end = s.find(del);

    while (end != -1) {
        splitted.push_back(s.substr(start, end - start));

        start = end + del.size();
        end = s.find(del, start);
    }

    splitted.push_back(s.substr(start, end - start));

    return splitted;
}

std::vector<Point> read_points(std::string filename) {
    std::vector<Point> points;
    std::ifstream File(filename);
    std::string line;

    while (std::getline(File, line)) {
        std::vector<std::string> point_strings = split(remove_characters(remove_characters(line, ')'), '('), " , ");

        for(int i = 0; i < point_strings.size(); i++) {
            std::vector<std::string> coordinate_strings = split(point_strings[i], ",");
            points.push_back(Point(std::stold(coordinate_strings[0]), std::stold(coordinate_strings[1])));
        }
    }

    File.close();

    return points;
}

void part1(list<Point> random_points) {
	clear_pixels();

//	list<Point> random_points = get_n_points_list(4);

//	write_and_plot("points.txt", random_points, 23);

	list<Point> closest_points = get_closest_points_bruteforce(random_points);

//	for (Point p : closest_points) {
//		plot_highlighted_point(p, 255, 0, 0);
//	}

//	write_to_ppm("output.ppm");

}

void part2(vector<Point> random_points) {
	clear_pixels();

//	vector<Point> random_points = get_n_points_vector(4);

	sort(random_points.begin(), random_points.end(), comparator);

//	for (Point p : random_points) {
//		cout << p.get_display_string() << endl;
//	}

	double min_dist = get_min_dist_recur(0, random_points.size(), random_points);

	cout << min_dist << endl;

}

int main() {
//	vector<Point> random_points_vector = read_points("points.txt");
	for (int i = 0; i < 50; i++) {
		vector<Point> random_points_vector = get_n_points_vector(60);
		list<Point> random_points_list;

		for (Point p : random_points_vector) {
			random_points_list.push_back(p);
		}

		cout << endl;

		part1(random_points_list);
		part2(random_points_vector);

		cout << endl;
	}

	return 0;
}
