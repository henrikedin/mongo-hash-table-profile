/* Measuring lookup times of unordered associative containers
 * without duplicate elements.
 *
 * Copyright 2013 Joaquin M Lopez Munoz.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 */

#include <algorithm>
#include <array>
#include <chrono>
#include <numeric>

//#include <boost/bind.hpp>
#include <iostream>
#include <random>
#include <vector>
#include "absl/strings/string_view.h"

#include <limits>       // std::numeric_limits
#include <algorithm>    // std::reverse
#include <thread>

#if defined(_MSC_VER)
#define FUNC_NOINLINE __declspec(noinline)
#else
#define FUNC_NOINLINE __attribute__ ((noinline))
#endif

namespace cppx {
	using std::numeric_limits;
	using std::reverse;

	typedef numeric_limits<long>    Long_info;
	int const long_digits = Long_info::max_digits10;
	int const long_bufsize = long_digits + 2;

	inline absl::string_view unsigned_to_decimal(unsigned long number, char* buffer)
	{
		char* start = buffer;
		if (number == 0)
		{
			*buffer++ = '0';
		}
		else
		{
			char* p_first = buffer;
			while (number != 0)
			{
				*buffer++ = '0' + number % 10;
				number /= 10;
			}
			reverse(p_first, buffer);
		}
		return absl::string_view(start, buffer - start);
	}

	inline auto decimal_from_unsigned(unsigned long number, char* buffer)
		-> absl::string_view
	{
		return unsigned_to_decimal(number, buffer);
	}

	inline absl::string_view to_decimal(long number, char* buffer)
	{
		if (number < 0)
		{
			buffer[0] = '-';
			return unsigned_to_decimal(-number, buffer + 1);
		}
		else
		{
			return unsigned_to_decimal(number, buffer);
		}
	}

	inline auto decimal_from(long number, char* buffer)
		-> absl::string_view
	{
		return to_decimal(number, buffer);
	}
}  // namespace cppx

using namespace std::chrono;
std::chrono::high_resolution_clock::time_point measure_start,measure_pause;

//template<typename F>
//double measure(F f)
//{
//  using namespace std::chrono;
//
//  static const int              num_trials=10;
//  static const milliseconds     min_time_per_trial(200);
//  std::array<double,num_trials> trials;
//
//  for(int i=0;i<num_trials;++i){
//    int                               runs=0;
//    high_resolution_clock::time_point t2;
//
//    measure_start=high_resolution_clock::now();
//    do{
//      f();
//      ++runs;
//      t2=high_resolution_clock::now();
//    }while(t2-measure_start<min_time_per_trial);
//    trials[i]=duration_cast<duration<double>>(t2 - measure_start).count()/runs;
//  }
//
//  std::sort(trials.begin(),trials.end());
//  return std::accumulate(
//    trials.begin()+2,trials.end()-2,0.0)/(trials.size()-4);
//}

template<typename F>
FUNC_NOINLINE double measure(F f, int n, int num_runs = 20)
{
	//static const int              num_trials = 10;
	//static const nanoseconds     min_time_per_trial(milliseconds(200));
	//std::array<double, num_trials> trials;

	//for (int i = 0; i<num_trials; ++i) {
	//	int                               runs = 0;
	//	high_resolution_clock::time_point t2;

	//	measure_start = high_resolution_clock::now();
	//	do {
	//		f();
	//		++runs;
	//		t2 = high_resolution_clock::now();
	//	} while (t2 - measure_start<min_time_per_trial);
	//	trials[i] = duration_cast<duration<double>>(t2 - measure_start).count() / runs;
	//}

	//std::sort(trials.begin(), trials.end());
	///*return duration_cast<nanoseconds>(duration<double>(
	//	std::accumulate(trials.begin() + 2, trials.end() - 2, 0.0) / (trials.size() - 4)));*/
	//return std::accumulate(trials.begin() + 2, trials.end() - 2, 0.0) / (trials.size() - 4);

	//int num_trials = (int)std::max(10ll, std::min(100000ll, 10000000ll - n_squared)); 
	int num_trials = std::min(2500000, n > 0 ? (int)ceil(25000000.0 / n) : 25000000);

	std::vector<double> runs;
	runs.resize(num_runs); 
	for (int r = 0; r < num_runs; ++r)
	{
		high_resolution_clock::time_point t2;
		measure_start = high_resolution_clock::now();
		for (int i = num_trials; i; --i)
		{
			f();
		}
		t2 = high_resolution_clock::now();
		runs[r] = duration_cast<duration<double>>(t2 - measure_start).count() / num_trials;
	}
	std::sort(runs.begin(), runs.end());
	return std::accumulate(runs.begin() + num_runs/5, runs.end() - num_runs / 5, 0.0) / (runs.size() - (num_runs / 5)*2);
}

void pause_timing()
{
  measure_pause=std::chrono::high_resolution_clock::now();
}

void resume_timing()
{
  measure_start+=std::chrono::high_resolution_clock::now()-measure_pause;
}

struct rand_seq
{
  rand_seq(unsigned int):gen(34862){}
  unsigned int operator()(){return dist(gen);}

private:
  std::uniform_int_distribution<unsigned int> dist;
  std::mt19937                                gen;
};

template<typename Container>
Container create(unsigned int n)
{
  Container s;
  rand_seq  rnd(n);
  while(n--)s.insert(rnd());
  return s;
}

template<typename Container>
struct uniform_distribution
{
	Container operator()(unsigned int n) const
	{
		return create<Container>(n);
	}
};

template<typename Container>
struct uniform_distribution_str
{
	Container operator()(unsigned int n) const
	{
		Container s;
		rand_seq  rnd(n);
		char str[64];
		while (n--)
		{
			s.emplace(cppx::decimal_from(rnd(), str));
		}
			
		return s;
	}
};

template<typename Container>
struct empty_container
{
	Container operator()(unsigned int n) const
	{
		Container s;
		return s;
	}
};

template<typename Container>
struct linear_sequence
{
	Container operator()(unsigned int n) const
	{
		Container s;
		for (unsigned int i = 0; i < n; ++i)
			s.insert(i);
		return s;
	}
};

template<typename Container>
struct scattered_successful_lookup
{
  typedef unsigned int result_type;

  unsigned int operator()(const Container & s,unsigned int n)const
  {
    unsigned int                                res=0;
    rand_seq                                    rnd(n);
    auto                                        end_=s.end();
    while(n--){
      if(s.find(rnd())!=end_)++res;
    }
    return res;
  }
};

template<typename Container>
struct scattered_successful_lookup_str
{
	typedef unsigned int result_type;

	unsigned int operator()(const Container & s, unsigned int n)const
	{
		unsigned int                                res = 0;
		rand_seq                                    rnd(n);
		char str[32];
		auto                                        end_ = s.end();
		while (n--) {
			if (s.find(cppx::decimal_from(rnd(), str)) != end_)++res;
		}
		return res;
	}
};

template<typename Container>
struct linear_successful_lookup
{
	typedef unsigned int result_type;

	unsigned int operator()(const Container & s, unsigned int n)const
	{
		unsigned int                                res = 0;
		auto                                        end_ = s.end();
		for (unsigned int i = 0; i < n; ++i) {
			if (s.find(i) != end_)++res;
		}
		return res;
	}
};

template<typename Container>
struct scattered_unsuccessful_lookup
{
  typedef unsigned int result_type;

  unsigned int operator()(const Container & s,unsigned int n)const
  {
    unsigned int                                res=0;
    std::uniform_int_distribution<unsigned int> dist;
    std::mt19937                                gen(76453);
    auto                                        end_=s.end();
    while(n--){
      if(s.find(dist(gen))!=end_)++res;
    }
    return res;
  }
};

template<typename Container>
void reserve(Container& s, unsigned int n)
{
	s.rehash((unsigned int)((double)n / s.max_load_factor()));
}

template<typename Container>
struct running_insertion
{
	typedef unsigned int result_type;

	unsigned int operator()(const Container &,unsigned int n)const
	{
		unsigned int res = 0;
		{
			Container s;
			rand_seq  rnd(n);
			while (n--)s.insert(rnd());
			res = s.size();
			pause_timing();
		}
		resume_timing();
		return res;
	}
};

template<typename Container>
struct linear_insertion
{
	typedef unsigned int result_type;

	unsigned int operator()(const Container &, unsigned int n)const
	{
		unsigned int res = 0;
		{
			Container s;
			for (unsigned int i=0; i < n; ++i)
				s.insert(i);
			res = s.size();
			pause_timing();
		}
		resume_timing();
		return res;
	}
};

template<typename Container>
struct norehash_running_insertion
{
	typedef unsigned int result_type;

	unsigned int operator()(unsigned int n)const
	{
		unsigned int res = 0;
		{
			Container s;
			rand_seq  rnd(n);
			reserve(s, n);
			while (n--)s.insert(rnd());
			res = s.size();
			pause_timing();
		}
		resume_timing();
		return res;
	}
};

template<typename K>
class noop_find
{
public:
	void* begin() { return this; }
	void* end() { return this; }
	const void* begin() const { return this; }
	const void* end() const { return this; }

	FUNC_NOINLINE void insert(const K& key) { }
	FUNC_NOINLINE void emplace(const K& key) { }
	FUNC_NOINLINE void* find(const K& key) { return this; }
	FUNC_NOINLINE const void* find(const K& key) const { return this; }
	std::size_t size() const { return 0; }

};


template<
  template<typename> class Tester,
  template<typename> class Creater,
  typename NoopContainer,
  typename Container1,typename Container2,typename Container3, typename Container4, typename Container5, typename Container6>
void test(
  const char* title,
  const char* name1,const char* name2,const char* name3, const char* name4, const char* name5, const char* name6)
{
  //unsigned int n0=1000,n1=3000000,dn=500;
  //double       fdn=1.025;

  unsigned int n0=500,n1=1000000;
  double       fdn=0.02;

  std::cout<<title<<":"<<std::endl;
  std::cout << ";"<<name1<<"; "<<name2<<"; "<<name3 << "; " << name4 << "; " << name5 << "; " << name6;
  std::cout << ";" << name1 << " (load_factor); " << name2 << " (load_factor); " << name3 << " (load_factor); " << name4 << " (load_factor); " << name5 << " (load_factor); " << name6 << " (load_factor); ";
  std::cout << std::endl;

  //for(unsigned int n=n0;n<=n1;n+=dn,dn=(unsigned int)(dn*fdn)){
  for (unsigned int n = n0; n <= n1; n += std::max(1u, (unsigned int)(n*fdn))) {
    double t = 0.0;
	double noop_t = 0.0;
	float lf[6];

	std::cout << n;
	{
		auto container_noop = Creater<NoopContainer>()(n);
		t = measure([&container_noop, n]() mutable
		{
			Tester<NoopContainer>()(container_noop, n);
		}, n, 40);
		noop_t = (t / n)*std::nano().den;
	}

	{
		auto container1 = Creater<Container1>()(n);
		t = measure([&container1, n]() mutable
		{
			Tester<Container1>()(container1, n);
		}, n);
		lf[0] = container1.load_factor();
		std::cout << ";" << ((t / n)*std::nano().den - noop_t);
	}
    
	{
		auto container2 = Creater<Container2>()(n);
		t = measure([&container2, n]() mutable
		{
			Tester<Container2>()(container2, n);
		}, n);
		lf[1] = container2.load_factor();
		std::cout << ";" << ((t / n)*std::nano().den - noop_t);
	}
	
	{
		auto container3 = Creater<Container3>()(n);
		t = measure([&container3, n]() mutable
		{
			Tester<Container3>()(container3, n);
		}, n);
		lf[2] = container3.load_factor();
		std::cout << ";" << ((t / n)*std::nano().den - noop_t);
	}

	{
		auto container4 = Creater<Container4>()(n);
		t = measure([&container4, n]() mutable
		{
			Tester<Container4>()(container4, n);
		}, n);
		lf[3] = container4.load_factor();
		std::cout << ";" << ((t / n)*std::nano().den - noop_t);
	}

	{
		auto container5 = Creater<Container5>()(n);
		t = measure([&container5, n]() mutable
		{
			Tester<Container5>()(container5, n);
		}, n);
		lf[4] = container5.load_factor();
		std::cout << ";" << ((t / n)*std::nano().den - noop_t);
	}

	{
		auto container6 = Creater<Container6>()(n);
		t = measure([&container6, n]() mutable
		{
			Tester<Container6>()(container6, n);
		}, n);
		lf[5] = container6.load_factor();
		std::cout << ";" << ((t / n)*std::nano().den - noop_t);
	}

	for (int i = 0; i < 6; ++i)
		std::cout << ";" << lf[i];

	std::cout << std::endl;
  }
}

  template<
	  template<typename> class Tester,
	  template<typename> class Creater,
	  typename NoopContainer,
	  typename Container1, typename Container2, typename Container3, typename Container4>
	  void test_4(
		  const char* title,
		  const char* name1, const char* name2, const char* name3, const char* name4)
  {
	  //unsigned int n0=1000,n1=3000000,dn=500;
	  //double       fdn=1.025;

	  unsigned int n0 = 500, n1 = 1000000;
	  double       fdn = 0.02;

	  std::cout << title << ":" << std::endl;
	  std::cout << ";" << name1 << "; " << name2 << "; " << name3 << "; " << name4;
	  std::cout << ";" << name1 << " (load_factor); " << name2 << " (load_factor); " << name3 << " (load_factor); " << name4 << " (load_factor); ";
	  std::cout << std::endl;

	  //for(unsigned int n=n0;n<=n1;n+=dn,dn=(unsigned int)(dn*fdn)){
	  for (unsigned int n = n0; n <= n1; n += std::max(1u, (unsigned int)(n*fdn))) {
		  double t;
		  double noop_t;
		  float lf[4];

		  std::cout << n;
		  {
			  auto container_noop = Creater<NoopContainer>()(n);
			  t = measure([&container_noop, n]() mutable
			  {
				  Tester<NoopContainer>()(container_noop, n);
			  }, n, 40);
			  noop_t = (t / n)*std::nano().den;
		  }

		  {
			  auto container1 = Creater<Container1>()(n);
			  t = measure([&container1, n]() mutable
			  {
				  Tester<Container1>()(container1, n);
			  }, n);
			  lf[0] = container1.load_factor();
			  std::cout << ";" << ((t / n)*std::nano().den - noop_t);
		  }

		  {
			  auto container2 = Creater<Container2>()(n);
			  t = measure([&container2, n]() mutable
			  {
				  Tester<Container2>()(container2, n);
			  }, n);
			  lf[1] = container2.load_factor();
			  std::cout << ";" << ((t / n)*std::nano().den - noop_t);
		  }

		  {
			  auto container3 = Creater<Container3>()(n);
			  t = measure([&container3, n]() mutable
			  {
				  Tester<Container3>()(container3, n);
			  }, n);
			  lf[2] = container3.load_factor();
			  std::cout << ";" << ((t / n)*std::nano().den - noop_t);
		  }

		  {
			  auto container4 = Creater<Container4>()(n);
			  t = measure([&container4, n]() mutable
			  {
				  Tester<Container4>()(container4, n);
			  }, n);
			  lf[3] = container4.load_factor();
			  std::cout << ";" << ((t / n)*std::nano().den - noop_t);
		  }


		  for (int i = 0; i < 4; ++i)
			  std::cout << ";" << lf[i];

		  std::cout << std::endl;
	  }
  }

  template<
	  template<typename> class Tester,
	  template<typename> class Creater,
	  typename NoopContainer,
	  typename Container1, typename Container2, typename Container3, typename Container4>
	  void test_4_empty(
		  const char* title,
		  const char* name1, const char* name2, const char* name3, const char* name4)
  {
	  //unsigned int n0=1000,n1=3000000,dn=500;
	  //double       fdn=1.025;

	  //unsigned int n0 = 500, n1 = 1000000;
	  //double       fdn = 0.02;

	  std::cout << title << ":" << std::endl;
	  std::cout << ";" << name1 << "; " << name2 << "; " << name3 << "; " << name4;
	  std::cout << ";" << name1 << " (load_factor); " << name2 << " (load_factor); " << name3 << " (load_factor); " << name4 << " (load_factor); ";
	  std::cout << std::endl;

	double t;
	double noop_t;
	unsigned int n = 1000;

	std::cout << 0;
	{
		auto container_noop = Creater<NoopContainer>()(0);
		t = measure([&container_noop, n]() mutable
		{
			Tester<NoopContainer>()(container_noop, n);
		}, n, 40);
		noop_t = (t/n)*std::nano().den;
	}

	{
		auto container1 = Creater<Container1>()(0);
		t = measure([&container1, n]() mutable
		{
			Tester<Container1>()(container1, n);
		}, n);
		std::cout << ";" << ((t / n)*std::nano().den - noop_t);
	}

	{
		auto container2 = Creater<Container2>()(0);
		t = measure([&container2, n]() mutable
		{
			Tester<Container2>()(container2, n);
		}, n);
		std::cout << ";" << ((t / n)*std::nano().den - noop_t);
	}

	{
		auto container3 = Creater<Container3>()(0);
		t = measure([&container3, n]() mutable
		{
			Tester<Container3>()(container3, n);
		}, n);
		std::cout << ";" << ((t / n)*std::nano().den - noop_t);
	}

	{
		auto container4 = Creater<Container4>()(0);
		t = measure([&container4, n]() mutable
		{
			Tester<Container4>()(container4, n);
		}, n);
		std::cout << ";" << ((t / n)*std::nano().den - noop_t);
	}

	std::cout << std::endl;
  }

//#include <cstdarg>
//#include <list>

struct result
{
	double time;
	float load_factor;
};

struct container_result
{
	std::string name;
	unsigned int n;
	std::vector<result> r;
};

//template<
//	template<typename> class Tester,
//	typename Container,
//	typename... Containers>
//	void test3(std::list<container_result>& results, va_list arglist)
//{
//	unsigned int n0 = 1000, n1 = 3000000, dn = 500;
//	double       fdn = 1.025;
//
//	container_result cr;
//	cr.name = va_arg(arglist, const char*);
//
//	auto& times = container_results.second;
//
//	for (unsigned int n = n0; n <= n1; n += dn, dn = (unsigned int)(dn*fdn)) {
//		double t;
//
//		auto container = create<Container>(n);
//		t = std::bind(Tester<Container>(), std::cref(container), n)();
//
//		result r;
//		r.time = (t / n)*10E6;
//		r.load_factor = container.load_factor();
//		times.push_back(r);
//	}
//	results.push_back(std::move(container_results));
//}
//
//template<
//	template<typename> class Tester,
//	typename... Containers>
//	void test2(
//		const char* title, ...)
//{
//	std::list<std::pair<std::string, std::vector<result>>> results;
//	va_list arglist;
//	va_start(arglist, title);
//	test3<Tester, Containers...>(results, arglist);
//	va_end(arglist);
//
//	std::cout << title << ":" << std::endl;
//	for (auto&& container : results)
//	{
//		std::cout << ";" << container.first;
//	}
//	std::cout << std::endl;
//
//	unsigned i = 0;
//	for (auto&& container : results)
//	{
//
//	}
//}

//template<
//	template<typename> class Tester,
//	typename Container1, typename Container2, typename Container3, typename Container4>
//	void test_insert(
//		const char* title,
//		const char* name1, const char* name2, const char* name3, const char* name4)
//{
//	unsigned int n0 = 1000, n1 = 3000000, dn = 500;
//	double       fdn = 1.025;
//
//	std::cout << title << ":" << std::endl;
//	std::cout << ";" << name1 << "; " << name2 << "; " << name3 << "; " << name4;
//	std::cout << std::endl;
//
//	for (unsigned int n = n0; n <= n1; n += dn, dn = (unsigned int)(dn*fdn)) {
//		double t;
//
//		std::cout << n;
//		{
//			t = measure([n]()
//			{
//				Tester<Container1>()(n);
//			});
//			std::cout << ";" << (t / n)*10E6;
//		}
//
//		{
//			t = measure([n]()
//			{
//				Tester<Container2>()(n);
//			});
//			std::cout << ";" << (t / n)*10E6;
//		}
//
//		{
//			t = measure([n]()
//			{
//				Tester<Container3>()(n);
//			});
//			std::cout << ";" << (t / n)*10E6;
//		}
//
//		{
//			t = measure([n]()
//			{
//				Tester<Container4>()(n);
//			});
//			std::cout << ";" << (t / n)*10E6;
//		}
//
//		std::cout << std::endl;
//	}
//}

//#include <boost/unordered_set.hpp>
//#include <boost/multi_index_container.hpp>
//#include <boost/multi_index/hashed_index.hpp>
//#include <boost/multi_index/identity.hpp>
#include <unordered_set>
#include <unordered_map>
#include "absl/container/flat_hash_set.h"
#include "absl/container/flat_hash_map.h"
#include "mongo/util/unordered_fast_key_table.h"
#include "bytell_hash_map.hpp"
#include "flat_hash_map.hpp"

#include "ska_modified/bytell_hash_map.hpp"
#include "murmurhash3/MurmurHash3.h"


template<typename T>
struct MongoTraits {
	static T hash(T a) {
		return a;
	}

	static bool equals(T a, T b) {
		return a == b;
	}

	static T toStorage(T s) {
		return s;
	}

	static T toLookup(T s) {
		return s;
	}

	class HashedKey {
	public:
		explicit HashedKey(T key = 0) : _key(key) {}

		HashedKey(T key, uint32_t hash) : _key(key) {
		}

		T key() const {
			return _key;
		}

		uint32_t hash() const {
			return _key;
		}

	private:
		T _key;
	};
};

struct AbslStringTraits {
	static uint32_t hash(absl::string_view a) {
		uint32_t hash;
		MurmurHash3_x86_32(a.data(), (int)a.size(), 0, &hash);
		return hash;
	}

	static bool equals(absl::string_view a, absl::string_view b) {
		return a == b;
	}

	static std::string toStorage(absl::string_view s) {
		return (std::string)s;
	}

	static absl::string_view toLookup(const std::string& s) {
		return absl::string_view(s);
	}

	class HashedKey {
	public:
		explicit HashedKey(absl::string_view key = "") : _key(key), _hash(AbslStringTraits::hash(_key)) {}

		HashedKey(absl::string_view key, uint32_t hash) : _key(key), _hash(hash) {
			// If you claim to know the hash, it better be correct.
		}

		const absl::string_view& key() const {
			return _key;
		}

		uint32_t hash() const {
			return _hash;
		}

	private:
		absl::string_view _key;
		uint32_t _hash;
	};
};

#include "absl/hash/internal/city.h"
struct AbslStringTraitsCity32 {
	static uint32_t hash(absl::string_view a) {
		return absl::hash_internal::CityHash64(a.data(), (int)a.size());
	}

	static bool equals(absl::string_view a, absl::string_view b) {
		return a == b;
	}

	static std::string toStorage(absl::string_view s) {
		return (std::string)s;
	}

	static absl::string_view toLookup(const std::string& s) {
		return absl::string_view(s);
	}

	class HashedKey {
	public:
		explicit HashedKey(absl::string_view key = "") : _key(key), _hash(AbslStringTraits::hash(_key)) {}

		HashedKey(absl::string_view key, uint32_t hash) : _key(key), _hash(hash) {
			// If you claim to know the hash, it better be correct.
		}

		const absl::string_view& key() const {
			return _key;
		}

		uint32_t hash() const {
			return _hash;
		}

	private:
		absl::string_view _key;
		uint32_t _hash;
	};
};

template <typename T>
struct NoopHasher {
	//using is_transparent = void;

	size_t operator()(T val) const {
		return val;
	}
};

struct payload128
{
	char p[128];
};

namespace std
{
	template<>
	struct hash<absl::string_view>
	{
		size_t operator()(const absl::string_view & a) const
		{
			uint32_t hash;
			MurmurHash3_x86_32(a.data(), (int)a.size(), 0, &hash);
			return hash;
		}
	};

}

struct absl_murmur_hash
{
	using is_transparent = void;

	size_t operator()(absl::string_view a) const
	{
		uint32_t hash;
		MurmurHash3_x86_32(a.data(), (int)a.size(), 0, &hash);
		return hash;
	}
};

int main()
{
  //typedef std::unordered_set<unsigned int>        container_t1;
  //typedef absl::flat_hash_set<unsigned int>       container_t2;
  //typedef absl::flat_hash_set<unsigned int, NoopHasher<unsigned int>>       container_t3;
  //typedef mongo::UnorderedFastKeyTable<unsigned int, unsigned int, bool, MongoTraits<unsigned int>> container_t4;
  //typedef ska::bytell_hash_set<unsigned int> container_t5;
  //typedef ska::flat_hash_set<unsigned int> container_t6;

  typedef absl::flat_hash_set<unsigned int>       container_t1;
  typedef mongo::UnorderedFastKeyTable<unsigned int, unsigned int, bool, MongoTraits<unsigned int>> container_t2;
  typedef ska::bytell_hash_set<unsigned int> container_t3;
  typedef absl::flat_hash_set<unsigned int, NoopHasher<unsigned int>>       container_t4;
  typedef std::unordered_set<unsigned int>        container_t5;
  typedef ska::flat_hash_set<unsigned int> container_t6;

  typedef absl::flat_hash_set<std::string>       container_t1_str;
  typedef mongo::UnorderedFastKeyTable<absl::string_view, std::string, bool, AbslStringTraits> container_t2_str;
  //typedef mongo::UnorderedFastKeyTable<absl::string_view, std::string, bool, AbslStringTraitsCity32> container_t3_str;
  typedef ska_modified::bytell_hash_set<absl::string_view, std::string>       container_t3_str;
  typedef absl::flat_hash_set<std::string, absl_murmur_hash>       container_t4_str;

  typedef std::unordered_map<unsigned int, payload128>        container_payload_t1;
  typedef absl::flat_hash_map<unsigned int, payload128>       container_payload_t2;
  typedef absl::flat_hash_map<unsigned int, NoopHasher<unsigned int>>       container_payload_t3;
  typedef mongo::UnorderedFastKeyTable<unsigned int, unsigned int, payload128, MongoTraits<unsigned int>> container_payload_t4;
  typedef ska::bytell_hash_map<unsigned int, payload128> container_payload_t5;
  typedef ska::flat_hash_map<unsigned int, payload128> container_payload_t6;

 test<
    scattered_successful_lookup,
	uniform_distribution,
	noop_find<unsigned int>,
    container_t1,
    container_t2,
    container_t3,
	container_t4,
	container_t5,
	container_t6>
  (
    "Scattered successful lookup",
	  "absl::flat_hash_set",
	  "mongo::UnorderedFastKeyTable",
	  "ska::bytell_hash_set",
	  "absl::flat_hash_set (noop hasher)",
	  "std::unordered_set",
	  "ska::flat_hash_set"
  );

  test<
    scattered_unsuccessful_lookup,
	uniform_distribution,
	noop_find<unsigned int>,
    container_t1,
    container_t2,
    container_t3,
    container_t4,
	container_t5,
	container_t6>
  (
    "Scattered unsuccessful lookup",
	  "absl::flat_hash_set",
	  "mongo::UnorderedFastKeyTable",
	  "ska::bytell_hash_set",
	  "absl::flat_hash_set (noop hasher)",
	  "std::unordered_set",
	  "ska::flat_hash_set"
  );

  test<
	  running_insertion,
	  empty_container,
	  noop_find<unsigned int>,
	  container_t1,
	  container_t2,
	  container_t3,
	  container_t4,
	  container_t5,
	  container_t6>
	(
	  "Runnning insertion",
		"absl::flat_hash_set",
		"mongo::UnorderedFastKeyTable",
		"ska::bytell_hash_set",
		"absl::flat_hash_set (noop hasher)",
		"std::unordered_set",
		"ska::flat_hash_set"
	);

  /*test<
	  scattered_unsuccessful_lookup,
	  linear_sequence,
	  noop_find<unsigned int>,
	  container_t1,
	  container_t2,
	  container_t3,
	  container_t4,
	  container_t5,
	  container_t6>
	(
		"Linear unsuccessful lookup",
		"absl::flat_hash_set",
		"mongo::UnorderedFastKeyTable",
		"ska::bytell_hash_set",
		"absl::flat_hash_set (noop hasher)",
		"std::unordered_set",
		"ska::flat_hash_set"
	);*/


  test_4<
	  scattered_successful_lookup_str,
	  uniform_distribution_str,
	  noop_find<absl::string_view>,
	  container_t1_str,
	  container_t2_str,
	  container_t3_str,
	  container_t4_str>
	  (
		  "Scattered successful lookup (string keys)",
		  "absl::flat_hash_set",
		  "mongo::UnorderedFastKeyTable",
		  "ska::bytell_hash_set",
		  "absl::flat_hash_set (murmur3)"
	  );

  test_4_empty<
	  scattered_successful_lookup_str,
	  uniform_distribution_str,
	  noop_find<absl::string_view>,
	  container_t1_str,
	  container_t2_str,
	  container_t3_str,
	  container_t4_str>
	  (
		  "Lookup in empty container (string keys)",
		  "absl::flat_hash_set",
		  "mongo::UnorderedFastKeyTable",
		  "ska::bytell_hash_set",
		  "absl::flat_hash_set (murmur3)"
	  );

  /*test<
	  scattered_unsuccessful_lookup,
	  uniform_distribution,
	  container_t1,
	  container_t2,
	  container_t3,
	  container_t4,
	  container_t5,
	  container_t6>
	  (
		  "Scattered unsuccessful lookup",
		  "std::unordered_set (MSVC 2015)",
		  "absl::flat_hash_set (default hasher)",
		  "absl::flat_hash_set (noop hasher)",
		  "mongo::UnorderedFastKeyTable",
		  "ska::bytell_hash_set",
		  "ska::flat_hash_set"
		  );

  test<
	  running_insertion,
	  empty_container,
	  container_t1,
	  container_t2,
	  container_t3,
	  container_t4,
	  container_t5,
	  container_t6>
	  (
		  "Runnning insertion",
		  "std::unordered_set (MSVC 2015)",
		  "absl::flat_hash_set (default hasher)",
		  "absl::flat_hash_set (noop hasher)",
		  "mongo::UnorderedFastKeyTable",
		  "ska::bytell_hash_set",
		  "ska::flat_hash_set"
		  );

  test<
	  scattered_unsuccessful_lookup,
	  linear_sequence,
	  container_t1,
	  container_t2,
	  container_t3,
	  container_t4,
	  container_t5,
	  container_t6>
	  (
		  "Linear unsuccessful lookup",
		  "std::unordered_set (MSVC 2015)",
		  "absl::flat_hash_set (default hasher)",
		  "absl::flat_hash_set (noop hasher)",
		  "mongo::UnorderedFastKeyTable",
		  "ska::bytell_hash_set",
		  "ska::flat_hash_set"
		  );*/

  /*test<
	  scattered_successful_lookup,
	  uniform_distribution,
	  container_payload_t1,
	  container_payload_t2,
	  container_payload_t3,
	  container_payload_t4,
	  container_payload_t5,
	  container_payload_t6>
	  (
		  "Scattered successful lookup (16 byte payload)",
		  "std::unordered_map (MSVC 2015)",
		  "absl::flat_hash_map (default hasher)",
		  "absl::flat_hash_map (noop hasher)",
		  "mongo::UnorderedFastKeyTable",
		  "ska::bytell_hash_map",
		  "ska::flat_hash_map"
		  );

  test<
	  scattered_unsuccessful_lookup,
	  uniform_distribution,
	  container_payload_t1,
	  container_payload_t2,
	  container_payload_t3,
	  container_payload_t4,
	  container_payload_t5,
	  container_payload_t6>
	  (
		  "Scattered unsuccessful lookup (16 byte payload)",
		  "std::unordered_map (MSVC 2015)",
		  "absl::flat_hash_map (default hasher)",
		  "absl::flat_hash_map (noop hasher)",
		  "mongo::UnorderedFastKeyTable",
		  "ska::bytell_hash_map",
		  "ska::flat_hash_map"
		  );

  test<
	  running_insertion,
	  empty_container,
	  container_payload_t1,
	  container_payload_t2,
	  container_payload_t3,
	  container_payload_t4,
	  container_payload_t5,
	  container_payload_t6>
	  (
		  "Runnning insertion (16 byte payload)",
		  "std::unordered_map (MSVC 2015)",
		  "absl::flat_hash_map (default hasher)",
		  "absl::flat_hash_map (noop hasher)",
		  "mongo::UnorderedFastKeyTable",
		  "ska::bytell_hash_map",
		  "ska::flat_hash_map"
		  );

  test<
	  scattered_unsuccessful_lookup,
	  linear_sequence,
	  container_payload_t1,
	  container_payload_t2,
	  container_payload_t3,
	  container_payload_t4,
	  container_payload_t5,
	  container_payload_t6>
	  (
		  "Linear unsuccessful lookup (16 byte payload)",
		  "std::unordered_map (MSVC 2015)",
		  "absl::flat_hash_map (default hasher)",
		  "absl::flat_hash_map (noop hasher)",
		  "mongo::UnorderedFastKeyTable",
		  "ska::bytell_hash_map",
		  "ska::flat_hash_map"
		  );*/
}
