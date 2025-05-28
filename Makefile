HTSLIB_DIR := ./external/htslib
CXX := $(if $(CXX),$(CXX),g++)
CXXFLAGS := -std=c++20 -Wall -Wextra
LDFLAGS := -L$(HTSLIB_DIR) -lhts -Wl,-rpath=$(abspath $(HTSLIB_DIR))
INCLUDES := -Iexternal/cxxopts/include -Iexternal/xtl/include -Iexternal/xtensor/include -I$(HTSLIB_DIR)

RELEASE_CXXFLAGS := -O3 -DNDEBUG -march=native
release: CXXFLAGS += $(RELEASE_CXXFLAGS)
release: all

SRC := src/transigner.cpp src/utils.cpp
OBJ := $(SRC:.cpp=.o)
TARGET := transigner

ifeq ($(shell $(CXX) -fopenmp -dM -E - < /dev/null | grep -q -i 'openmp'; echo $$?),0)
    OPENMP := -fopenmp
    CXXFLAGS += $(OPENMP)
endif

all: $(HTSLIB_DIR)/libhts.a $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(HTSLIB_DIR)/libhts.a:
	cd $(HTSLIB_DIR) && autoreconf -i
	cd $(HTSLIB_DIR) && ./configure
	cd $(HTSLIB_DIR) && $(MAKE)

clean:
	rm -f $(OBJ) $(TARGET)
	$(MAKE) -C $(HTSLIB_DIR) clean