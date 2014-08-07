package varbin;

public enum Genotype {
	WILDTYPE, // GT = 0/0
	HET, // GT = 0/1
	HOM, // GT = 1/1
	NOCALL // GT = ./. (not enough read info to make a genotype call)

}
