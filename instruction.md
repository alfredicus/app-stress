# Instructions

## Installation of the tools
Several tools have to be installed:

### 1. nodejs

### 2. npm

### 3. yarn

### 4. Git

### 5. Visual Studio Code

## GitHub
Explanation...

## npmjs
Explanation...

## Organization of the codes 

## Debugging with jest
A la racine de la lib stress:
```sh
node --inspect-brk node_modules/.bin/jest --runInBand src/tests/stress-tensor.test.ts
```

Puis dans Chrome, dans la barre d’adresse, taper:
```text
chrome://inspect
```
Et choisir le lien

## Extending the type of curve
The interface is:
```ts
export interface GenericCurve {
    generate(theta: number, phi: number): string
}
```

