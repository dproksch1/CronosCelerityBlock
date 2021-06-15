#include "buildinfo.H"

#ifndef CRONOS_GIT_VERSION
# define CRONOS_GIT_VERSION "N/A"
#endif

#ifndef CRONOS_GIT_COMMIT
# define CRONOS_GIT_COMMIT "N/A"
#endif

#ifndef APP_GIT_REPO
# define APP_GIT_REPO "N/A"
#endif

#ifndef APP_GIT_VERSION
# define APP_GIT_VERSION "N/A"
#endif

#ifndef APP_GIT_COMMIT
# define APP_GIT_COMMIT "N/A"
#endif

#ifndef BUILD_DATE
# define BUILD_DATE "N/A"
#endif

#ifndef CRONOS_ROOT_DIR
# define CRONOS_ROOT_DIR "N/A"
#endif


std::string CronosGitVersion() {
	return CRONOS_GIT_VERSION;
}

std::string CronosGitCommit() {
	return CRONOS_GIT_COMMIT;
}

std::string ApplicationGitRepo() {
	return APP_GIT_REPO;
}

std::string ApplicationGitVersion() {
	return APP_GIT_VERSION;
}

std::string ApplicationGitCommit() {
	return APP_GIT_COMMIT;
}

std::string BuildDate() {
	return BUILD_DATE;
}

std::string CronosRootDir() {
	return CRONOS_ROOT_DIR;
}
