#ifndef BUILDINFO_H
#define BUILDINFO_H

#include <string>

std::string CronosGitVersion();
std::string CronosGitCommit();
std::string ApplicationGitRepo();
std::string ApplicationGitVersion();
std::string ApplicationGitCommit();
std::string BuildDate();
std::string CronosRootDir();

#endif
