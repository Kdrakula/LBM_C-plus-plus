CXX = c++
CXXFLAGS = -std=c++17 -O3 -Wall

TARGET = lbm_simulation
SRC = main.cpp
OUT_DIR = simulation_output
VIDEO = simulation.mp4
FRAMERATE = 10

all: video

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

build: $(TARGET)

run: build
	rm -rf $(OUT_DIR)
	mkdir -p $(OUT_DIR)
	./$(TARGET)

video: run
	ffmpeg -y -framerate $(FRAMERATE) -pattern_type glob -i '$(OUT_DIR)/*.ppm' -c:v libx264 -pix_fmt yuv420p $(VIDEO)
clean:
	rm -f $(TARGET) $(VIDEO)
	rm -rf $(OUT_DIR)

.PHONY: all build run video clean