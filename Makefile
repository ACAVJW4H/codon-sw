.PHONY: debug release compile setup-cmake clean style doc

BUILD = _build

release: CMAKE_BUILD_TYPE = Release
release: BUILD_DIR = $(BUILD)/$(CMAKE_BUILD_TYPE)
release: compile

debug: CMAKE_BUILD_TYPE = Debug
debug: BUILD_DIR = $(BUILD)/$(CMAKE_BUILD_TYPE)
debug: compile

compile: setup-cmake
	+make -C$(BUILD_DIR)

setup-cmake:
	mkdir -p $(BUILD_DIR)
	cd $(BUILD_DIR) && cmake -DCMAKE_BUILD_TYPE=$(CMAKE_BUILD_TYPE) ../..

clean:
	rm -rf $(BUILD) html latex

doc:
	doxygen

style:
	astyle  -A3 \
	        --pad-oper \
	        --unpad-paren \
	        --keep-one-line-blocks \
	        --keep-one-line-statements \
	        --suffix=none \
	        --formatted \
	        --lineend=linux \
					--align-pointer=type \
	        `find src -regextype posix-extended -regex ".*\.(cc|h|hpp|cpp)$$"`
