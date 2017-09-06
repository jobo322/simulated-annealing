import Matrix from "ml-matrix";

function simulatedAnnealing(inputs) {
    var {
        objetiveFunction,
        guess,
        lowerBound,
        upperBound,
        iterationsNumber,
        quenchingFactor,
        toleranceValue
    } = inputs;

    let boundsDiff = upperBound.clone().subtract(lowerBound);
    let x = guess;
    let fx = objetiveFunction(x);
    let xi = x;
    let fi = fx;
    let initialErrorFunctionValue = fi;
    let acceptableError = initialErrorFunctionValue * 0.1;
    let j = 0;
    for (let k = 0; k < iterationsNumber; k++) {
        let ti =   Math.pow(k / iterationsNumber, quenchingFactor);
        let mu = Math.pow(10, ti*100 );
        let dx = muInv(Matrix.rand(x.rows, x.columns).multiply(2).add(-1), mu).multiply(boundsDiff);
        let x1 = x.clone().add(dx);
        for (let i = 0; i < x1.columns; i ++) {
            x1[0][i] = (x1[0][i] < lowerBound[0][i] ? lowerBound[0][i] : 0) +
            (lowerBound[0][i] <= x1[0][i] && x1[0][i] <= upperBound[0][i] ? x1[0][i] : 0) +
            (upperBound[0][i] < x1[0][i] ? upperBound[0][i] : 0);
        }
        let fx1 = objetiveFunction(x1);
        let df = fx1 - fx;
        if(df < 0 || Math.random < Math.random < (Math.exp((-ti * df) / (Math.abs(fx) + 8E-12) / toleranceValue))) {
            for (var i = 0; i < x.columns; i++) {
                x[0][i] = x1[0][i];
            }
            fx = fx1;
        }
        if(fx < fi) {
            fi = fx1;
            for (let i = 0; i < x.columns; i++) {
                xi[0][i] = x[0][i];
            }
        }
        if(fi < acceptableError) {
            k = iterationsNumber;
        }
    }
    return [xi, fi];
}

function muInv(y, mu){
    let x = Matrix.zeros(y.rows, y.columns);
    for (let i = 0; i < y.columns; i++) {
        x[0][i] = (((Math.pow(1 + mu, Math.abs(y[0][i])) - 1 ) / mu)) * Math.sign(y[0][i]);
    }
    return x;
}

module.exports = simulatedAnnealing;
