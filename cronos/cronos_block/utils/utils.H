#ifndef CRONOS_UTILS_H
#define CRONOS_UTILS_H = 1

#define __xstr(a) __insieme_utils_str(a)
#define __insieme_utils_str(a) #a

#ifdef __GNUC__
#define __insieme_unused __attribute__((unused))
#else
#define __insieme_unused
#endif

namespace insieme {
	namespace utils {

		// allows to provide a lambda that adds printed information when an assertion occurs
		//void setAssertExtraInfoPrinter(std::function<void(void)> printer);
		// clears the extra information printer lambda
		//void clearAssertExtraInfoPrinter();

		namespace detail {

			//extern thread_local std::function<void(void)> extraAssertInformationPrinter;

			struct LazyAssertion {
				bool value;
				LazyAssertion(bool value) : value(value) {}
				~LazyAssertion() /*{ if(!value) __debugbreak(); } //*/= default;
				operator bool() const {
					return !value;
				}
			};

		} // end namespace detail
	} // end namespace utils
} // end namespace insieme

#define assert_fail()                                                                                                                                          \
	if(__insieme_unused auto __x = insieme::utils::detail::LazyAssertion(false)) std::cerr << "\nAssertion failed in " __FILE__ ":" __xstr(__LINE__) " - "


namespace gpu {

	enum {

		FaceWest = 0,
		FaceEast,
		FaceSouth,
		FaceNorth,
		FaceBottom,
		FaceTop,

		FaceMax,

	};

}

#endif