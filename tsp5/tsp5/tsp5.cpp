#include <iostream>
#include <fstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <chrono>
#include <iomanip>
#include <random>

struct Point {
    double x;
    double y;
};

double distance(const Point& p1, const Point& p2) {
    return std::sqrt(std::pow(p1.x - p2.x, 2) + std::pow(p1.y - p2.y, 2));
}

int fitness(const std::vector<int>& tour, const std::vector<Point>& points) {
    int total_distance = 0;
    for (size_t i = 0; i < tour.size() - 1; ++i) {
        total_distance += distance(points[tour[i]], points[tour[i + 1]]);
    }
    total_distance += distance(points[tour.back()], points[tour.front()]); // Return to the starting point
    return total_distance;
}

std::vector<int> random_neighbor(const std::vector<int>& tour, std::mt19937& gen) {
    std::vector<int> neighbor = tour;
    std::uniform_int_distribution<> dis(1, neighbor.size() - 2); // Exclude the first and last cities
    int i = dis(gen);
    int j = dis(gen);
    std::reverse(neighbor.begin() + std::min(i, j), neighbor.begin() + std::max(i, j));
    return neighbor;
}

std::vector<int> nearest_neighbor_tour(const std::vector<Point>& points) {
    std::vector<int> tour(points.size());
    std::iota(tour.begin(), tour.end(), 0);

    for (size_t i = 0; i < points.size() - 1; ++i) {
        double min_dist = std::numeric_limits<double>::max();
        size_t nearest_idx = i + 1;
        for (size_t j = i + 1; j < points.size(); ++j) {
            double dist = distance(points[tour[i]], points[tour[j]]);
            if (dist < min_dist) {
                min_dist = dist;
                nearest_idx = j;
            }
        }
        std::swap(tour[i + 1], tour[nearest_idx]);
    }

    return tour;
}

std::vector<int> simulated_annealing(const std::vector<Point>& points, int max_iterations, double initial_temperature, double cooling_rate, double& optimal_distance) {
    std::vector<int> current_tour = nearest_neighbor_tour(points);
    std::vector<int> best_tour = current_tour;
    double current_distance = fitness(current_tour, points);
    double best_distance = current_distance;

    std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<> dis(0, 1);

    for (int iter = 0; iter < max_iterations; ++iter) {
        std::vector<int> neighbor = random_neighbor(current_tour, gen);
        double neighbor_distance = fitness(neighbor, points);

        double temperature = initial_temperature * std::pow(cooling_rate, iter);
        double acceptance_probability = std::exp((current_distance - neighbor_distance) / temperature);

        double random_value = dis(gen);

        if (neighbor_distance < current_distance || random_value < acceptance_probability) {
            current_tour = neighbor;
            current_distance = neighbor_distance;
        }

        if (current_distance < best_distance) {
            best_tour = current_tour;
            best_distance = current_distance;
        }
    }

    optimal_distance = best_distance;
    return best_tour;
}

std::vector<Point> read_points_from_file(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<Point> points;
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return points;
    }
    int num_points;
    file >> num_points;
    points.resize(num_points);
    for (int i = 0; i < num_points; ++i) {
        file >> points[i].x >> points[i].y;
    }
    file.close();
    return points;
}

int main() {
    auto start = std::chrono::high_resolution_clock::now();
    std::string filename = "c:/Users/LENOVO/Desktop/Data/tsp_85900_1";
    std::vector<Point> points = read_points_from_file(filename);
    if (points.empty()) {
        std::cerr << "Error: No points read from file." << std::endl;
        return 1;
    }

    int max_iterations = 5000;
    double initial_temperature = 1000.0;
    double cooling_rate = 0.999;

    double optimal_distance;
    std::vector<int> best_tour = simulated_annealing(points, max_iterations, initial_temperature, cooling_rate, optimal_distance);

    std::cout << "Best tour: ";
    for (int city : best_tour) {
        std::cout << city << " ";
    }
    std::cout << std::endl;

    std::cout << "Optimal distance: " << std::fixed << std::setprecision(2) << optimal_distance << std::endl;
    auto end = std::chrono::high_resolution_clock::now();
    auto duration_seconds = std::chrono::duration<double>(end - start).count();
    std::cout << "Program running time: " << std::fixed << std::setprecision(2) << duration_seconds << " seconds" << std::endl;

    return 0;
}
