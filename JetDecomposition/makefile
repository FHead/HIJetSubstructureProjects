default: CompileAll

CompileAll: bin/BasicJetImageAnalysis bin/TruncationTest

RunAll: RunBasicJetImageAnalysis RunTruncationTest

bin/BasicJetImageAnalysis: source/BasicJetImageAnalysis.cpp
	g++ source/BasicJetImageAnalysis.cpp -o ./bin/BasicJetImageAnalysis \
		-I./include \
		`root-config --glibs --cflags` \
		`fastjet-config --cxxflags --libs`

RunBasicJetImageAnalysis: bin/BasicJetImageAnalysis
	mkdir -p output
	./bin/BasicJetImageAnalysis --input Samples/JetTreeSample.root \
		--tree DiscretizedJetTree \
		--output output/JetImages \
		--minjetpt 30

bin/FTHelper.o: source/FTHelper.cpp include/FTHelper.h
	g++ source/FTHelper.cpp -c -o bin/FTHelper.o -Iinclude -std=c++17

bin/TruncationTest: source/TruncationTest.cpp bin/FTHelper.o
	g++ source/TruncationTest.cpp -o ./bin/TruncationTest \
		-Iinclude \
		`root-config --glibs --cflags` \
		bin/FTHelper.o

RunTruncationTest: bin/TruncationTest
	mkdir -p output
	./bin/TruncationTest --Input Samples/JetTree_0.root --Output output/FourierTransformTruncation.pdf \
		--Percentage 0.1 --DoIndividual true --MinS 100
