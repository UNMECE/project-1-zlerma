// Zaid Orozco-Lerma
// Project 1

#include <iostream>
#include <iomanip> // This is something I to help with the formatting of the output
#include <math> // I don't think I ended up needing this but I gotta submit and I don't want to mess it up now

// Here is the struct, as was specified in the assignment
struct _capacitor
{
	double *time;
	double *voltage;
	double *current;
	double C;
};

typedef struct _capacitor Capacitor;

// A whole bunch of constants
const double delta_t = 1e-10;
const double final_time = 5e-6;
// This is to adjust num_timesteps to include the final time point
const int num_timesteps = static_cast<int>(final_time / delta_t) + 1; // I forgot why I did +1, but it's for an error reason
const double R = 1000; // 1 kΩ
const double C_value = 100e-12; // 100 pF
const double I_const = 1e-2; // 10 mA
const double V0 = 10.0; // Constant voltage case

// Function prototypes
void constant_current_simulation(Capacitor& capacitor);
void constant_voltage_simulation(Capacitor& capacitor);
void allocate_memory(Capacitor& capacitor);
void free_memory(Capacitor& capacitor);

int main()
{
	// Create a Capacitor object
	Capacitor capacitor;
	capacitor.C = C_value; // Set the capacitance to the value we defined earlier

	// Allocate memory for the arrays in the capacitor struct
	allocate_memory(capacitor);

	// Fill the time array with time steps
	for (int i = 0; i < num_timesteps; ++i) {
		capacitor.time[i] = i * delta_t; // Each time point is just i times delta_t
	}

	// Simulate the constant current charging scenario
	std::cout << "Simulating constant current charging..." << std::endl;
	constant_current_simulation(capacitor);

	// Simulate the constant voltage charging scenario
	std::cout << "\nSimulating constant voltage charging..." << std::endl;
	constant_voltage_simulation(capacitor);

	// Free up the memory we allocated
	free_memory(capacitor);

	// Successful! If we get here
	return 0;
}

// Function to simulate the constant current charging scenario
void constant_current_simulation(Capacitor& capacitor)
{
	// Initial conditions
	capacitor.voltage[0] = 0.0;	// Start with zero voltage across the capacitor
	capacitor.current[0] = I_const; // Current is constant at I_const

	// Loop over each time step to calculate voltage (a LOT of times)
	for (int i = 1; i < num_timesteps; ++i) {
		// Current remains constant in this scenario!
		capacitor.current[i] = I_const;

		// Calculate the voltage at the next time step using the finite-difference method
		// V(t+1) = V(t) + (I(t) * delta_t) / C
		capacitor.voltage[i] = capacitor.voltage[i - 1] + (capacitor.current[i - 1] * delta_t) / capacitor.C;

		// Output every 200 time steps
		if (i % 200 == 0) {
			std::cout << std::fixed << std::setprecision(6); // Set precision for better readability
			std::cout << "Time: " << std::setw(10) << capacitor.time[i]
					  << " s, Current: " << std::setw(10) << capacitor.current[i]
					  << " A, Voltage: " << std::setw(10) << capacitor.voltage[i] << " V" << std::endl;
		}
	}
}

// Function to simulate the constant voltage charging scenario
void constant_voltage_simulation(Capacitor& capacitor)
{
	// Initial conditions
	capacitor.current[0] = V0 / R; // I(0) = V0 / R, based on Ohm's law
	capacitor.voltage[0] = 0.0;	 // Start with zero voltage across the capacitor

	// Loop over each time step to calculate current and voltage:
	for (int i = 1; i < num_timesteps; ++i) {
		// Calculate the current at the next time step:
		// I(t+1) = I(t) - (I(t) / (R * C)) * delta_t
		capacitor.current[i] = capacitor.current[i - 1] - (capacitor.current[i - 1] / (R * capacitor.C)) * delta_t;

		// Calculate the voltage at the next time step:
		// V(t+1) = V(t) + (I(t) * delta_t) / C
		capacitor.voltage[i] = capacitor.voltage[i - 1] + (capacitor.current[i - 1] * delta_t) / capacitor.C;

		// Output every 200 time steps
		if (i % 200 == 0) {
			std::cout << std::fixed << std::setprecision(6); // Again, setting precision
			std::cout << "Time: " << std::setw(10) << capacitor.time[i]
					  << " s, Current: " << std::setw(10) << capacitor.current[i]
					  << " A, Voltage: " << std::setw(10) << capacitor.voltage[i] << " V" << std::endl;
		}
	}
}

// We have to allocate memory for the arrays
void allocate_memory(Capacitor& capacitor)
{
	// Allocate memory for time, voltage, and current arrays
	capacitor.time = new double[num_timesteps];
	capacitor.voltage = new double[num_timesteps];
	capacitor.current = new double[num_timesteps];

	// Check if memory allocation was successful
	if (!capacitor.time || !capacitor.voltage || !capacitor.current) {
		std::cerr << "Memory allocation failed!!! :(" << std::endl;
		exit(1); // Exit the program if allocation fails
	}

	// Initialize the arrays to zero, just to be safe
	for (int i = 0; i < num_timesteps; ++i) {
		capacitor.time[i] = 0.0;
		capacitor.voltage[i] = 0.0;
		capacitor.current[i] = 0.0;
	}
}

// Function to free the dynamically allocated memory
void free_memory(Capacitor& capacitor)
{
	// Delete the arrays to free up memory
	delete[] capacitor.time;
	delete[] capacitor.voltage;
	delete[] capacitor.current;
}
