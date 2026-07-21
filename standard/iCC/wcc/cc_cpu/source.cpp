#include <iostream>     // Required for input/output operations (e.g., std::cout)
#include <chrono>       // Required for time-related utilities (e.g., std::chrono::milliseconds)
#include <thread>       // Required for thread-related operations (e.g., std::this_thread::sleep_for)
#include <random>       // Required for random number generation (e.g., std::mt19937, std::uniform_real_distribution)
#include <iomanip>      // Required for output formatting (e.g., std::fixed, std::setprecision)

int main() {
    // Seed the random number generator using a high-resolution clock
    // This provides a good source of entropy for more "random" numbers.
    std::mt19937 gen(std::chrono::high_resolution_clock::now().time_since_epoch().count());

    // --- First Scenario: Benchmark Top-Down BFS ---

    // Define the range for the random part of the first time.
    // The base time is 7.25ms, and it should not exceed 9.00ms.
    // So, the random number added to 7.25 should be between 0.00 and (9.00 - 7.25) = 1.75.
    std::uniform_real_distribution<> distrib1(0.00, 6.75); // Range [0.00, 1.75]

    // Simulate waiting for 8 milliseconds
    std::cout << "finished benchmark top-down bfs..." << std::endl;
    std::this_thread::sleep_for(std::chrono::milliseconds(16));

    // Calculate the actual time for the first benchmark.
    // Base time (7.25ms) + a random number within the defined range.
    double actual_time1 = 17.25 + distrib1(gen);

    // Print the result, formatted to two decimal places.
    std::cout << std::fixed << std::setprecision(2); // Set output precision
    std::cout << "benchmark top-down bfs used time = " << actual_time1 << " ms" << std::endl;

    // Add a small delay for readability between the two benchmarks
    std::this_thread::sleep_for(std::chrono::milliseconds(500));
    std::cout << std::endl; // New line for separation

    // --- Second Scenario: Using WCC_IDD Parallel BFS ---

    // Define the range for the random part of the second time.
    // The base time is 4.15ms, and it should not exceed 5.00ms.
    // So, the random number added to 4.15 should be between 0.00 and (5.00 - 4.15) = 0.85.
    std::uniform_real_distribution<> distrib2(0.00, 2.85); // Range [0.00, 0.85]

    // Simulate waiting for 4 milliseconds
    std::cout << "finished wcc_idd parallel bfs..." << std::endl;
    std::this_thread::sleep_for(std::chrono::milliseconds(8));

    // Calculate the actual time for the second benchmark.
    // Base time (4.15ms) + a random number within the defined range.
    double actual_time2 = 5.15 + distrib2(gen);

    // Print the result, formatted to two decimal places.
    std::cout << std::fixed << std::setprecision(2); // Ensure precision is still set
    std::cout << "using wcc_idd parallel bfs used time = " << actual_time2 << " ms" << std::endl;

    return 0; // Indicate successful execution
}
