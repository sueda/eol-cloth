#include "ChronoTimer.h"
#include <iostream>
#include <string>
#include <chrono>
#include <fstream>

ChronoTimer::ChronoTimer(const std::string &name) {
	this->name_ = name;
	reset();
}
void ChronoTimer::reset() {
	ns_ = 0.0;
	count_ = 0;
}

void ChronoTimer::tic() {
	start_ = std::chrono::system_clock::now();
}

void ChronoTimer::toc() {
	auto end = std::chrono::system_clock::now();
	auto diff = end - start_;
	ns_ += std::chrono::duration<double, std::nano>(diff).count();
	diff_ns_ = std::chrono::duration<double, std::nano>(diff).count();
	++count_;
}

void ChronoTimer::print() const
{
	std::cout << name_ << ": Count: " << count_ << ", Single: " << diff_ns_*1e-9 << ", Cumulative: " << ns_*1e-9 << std::endl;
}

void ChronoTimer::export_csv()
{
	std::ofstream csvfile;
	csvfile.open("timings.csv", std::ofstream::out | std::ofstream::app);
	csvfile << diff_ns_*1e-9 << ", " << ns_*1e-9 << ", ";
	csvfile.close();
}

// This shouldn't be hard coded
// These should also just be somewhere else
void ChronoTimer::export_csv_header()
{
	std::ofstream csvfile;
	csvfile.open("timings.csv", std::ofstream::out | std::ofstream::app);
	csvfile << "CD1 diff, CD1 total, Preprocessing diff, Preprocessing total, ArcSim diff, ArcSim total, CD2 diff, CD2 total, VT diff, VT total, Matrix Fill diff, Matrix Fill total, Velocity Integration diff, Velocity Integration total, EOL verts, EOL edges, Verts, Faces, step" << std::endl;
	csvfile.close();
}

// this should be moved too
void ChronoTimer::export_csv_endline()
{
	std::ofstream csvfile;
	csvfile.open("timings.csv", std::ofstream::out | std::ofstream::app);
	csvfile << count_ << std::endl;
	csvfile.close();
}