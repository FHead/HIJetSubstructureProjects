default: CompileAll

CompileAll: Bin/BasicJetImageAnalysis

RunAll: RunBasicJetImageAnalysis

Bin/BasicJetImageAnalysis: Source/BasicJetImageAnalysis.cpp
	g++ Source/BasicJetImageAnalysis.cpp -o ./Bin/BasicJetImageAnalysis \
		-I./Include \
		`root-config --glibs --cflags` \
		`fastjet-config --cxxflags --libs`

RunBasicJetImageAnalysis: Bin/BasicJetImageAnalysis
	./Bin/BasicJetImageAnalysis --input Samples/JetTreeSample.root \
		--tree DiscretizedJetTree \
		--output JetImages \
		--minjetpt 30
