#pragma once
#ifndef __ChronoTimer__
#define __ChronoTimer__

#include <iostream>
#include <chrono>


class ChronoTimer
{
private:
	std::string name_;
	double ns_;
	double diff_ns_;
	int count_;
	std::chrono::time_point<std::chrono::system_clock> start_;

public:
	ChronoTimer(const std::string &name);
	void reset();
	void tic();
	void toc();
	void print() const;
	void export_csv();
	void export_csv_header();
	void export_csv_endline();
};

#endif