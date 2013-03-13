package seqSimulator.evolution.proteinstruct;

import java.util.ArrayList;
import java.util.List;

public class StructureEnv {

	List<double [][]> structEnv;
	List<double [][]> logStructEnv;
	
	// first row represents row marginal, second row represents column marginal
	List<double [][]> marginalProbMatrices;
	List<double [][]> marginalLogProbMatrices;
	
	public StructureEnv(){
		structEnv = new ArrayList<double[][]>();
	}
	
	public void addEnv(double[][] newEnv) {
		structEnv.add(newEnv);
	}
	
	public double getLogProb(int structEnvNumber, int firstCodonType,
			int secondCodonType) {
		return logStructEnv.get(structEnvNumber)[firstCodonType][secondCodonType];
	}

	public void calculateLogStructEnv(){
		this.getMarginalProbMatrices();
		this.calculateMarginalLogProbMatrices();
		
		logStructEnv  = new ArrayList<double [][]>();
		// take each env and pre-compute log value for each matrix
		for (int envNum = 0; envNum < structEnv.size(); envNum++) {
			logStructEnv.add(getLogMatrixAdjustedByMarginal(structEnv.get(envNum), marginalLogProbMatrices.get(envNum)));
		}
	}
	
	// getter
	void getMarginalProbMatrices () {
		marginalProbMatrices = new ArrayList<double [][]>();
		for (int envNum = 0; envNum < structEnv.size(); envNum++) {
			marginalProbMatrices.add(getMarginalProbMatrix(structEnv.get(envNum)));
		}
	}
	
	public double [][] getMarginalProbMatrix(double[][] structEnvMatrix){
		int structEnvMatrixDim = structEnvMatrix.length;
		// create [2][61] matrix
		double [][] marginalMatrix = new double [2][structEnvMatrixDim];
		// calculate row & col marginal
		for (int rowNr = 0 ; rowNr < structEnvMatrixDim; rowNr++) {
			marginalMatrix[0][rowNr] = getRowSum(rowNr, structEnvMatrix);
			int colNr = rowNr;
			marginalMatrix[1][colNr] = getColSum(colNr, structEnvMatrix);
		}
		// return result
		return marginalMatrix;
	}
	
	double getRowSum (int rowNr, double[][] structEnvMatrix) {
		double rowSum = 0;
		for (int colNr = 0; colNr < structEnvMatrix[0].length; colNr++) {
			rowSum += structEnvMatrix[rowNr][colNr];
		}
		return rowSum;
	}
	
	double getColSum (int colNr, double[][] structEnvMatrix) {
		double colSum = 0;
		for (int rowNr = 0; rowNr < structEnvMatrix.length; rowNr++) {
			colSum += structEnvMatrix[rowNr][colNr];
		}
		return colSum;
	}
	
	void calculateMarginalLogProbMatrices(){
		marginalLogProbMatrices  = new ArrayList<double [][]>();
		
		for (int envNum = 0; envNum < marginalProbMatrices.size(); envNum++) {
			marginalLogProbMatrices.add(getLogMatrix(marginalProbMatrices.get(envNum)));
		}
	}
	
	double [][] getLogMatrix(double [][] matrix){
		int rowDim = matrix.length;
		int colDim = matrix[0].length;
		double [][] logMatrix = new double[rowDim][colDim];
		
		for(int i = 0; i < rowDim; i++){
			for (int j = 0 ;j< colDim; j++){
				logMatrix[i][j] = Math.log(matrix[i][j]);
			}
		}
		
		return logMatrix;
	}
	
	double [][] getLogMatrixAdjustedByMarginal(double [][] matrix, double[][] marginalLogMatrix){
		int rowDim = matrix.length;
		int colDim = matrix[0].length;
		double [][] logMatrix = new double[rowDim][colDim];
		
		for(int i = 0; i < rowDim; i++){
			for (int j = 0 ;j< colDim; j++){
				logMatrix[i][j] = Math.log(matrix[i][j]) - marginalLogMatrix[0][i] - marginalLogMatrix[1][j];
			}
		}
		
		return logMatrix;
	}	
	
	// for debugging
	public double getEnvEntry(int envNr, int i, int j){
		return structEnv.get(envNr)[i][j];
	}
	
	public double getMatrixSum(int envNr){
		double sum = 0;
		double[][] matrix = structEnv.get(envNr);
		for(int i = 0; i < 61; i++){
			for (int j = 0; j < 61; j++){
				sum += matrix[i][j];
			}
		}
		return sum;
	}
}
