#pragma once
#ifndef __parseParams__
#define __parseParams__

#include "genSet.h"
#include "Scene.h"
#include <string>
#include <memory>

void load_genset(const std::shared_ptr<genSet> genset, const std::string &JSON_FILE);

void load_simset(std::shared_ptr<Scene> scene, const std::string &JSON_FILE);

#endif
