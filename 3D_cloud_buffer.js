import * as THREE from 'https://unpkg.com/three@0.126.1/build/three.module.js';

import { OrbitControls } from 'https://unpkg.com/three@0.126.1/examples/jsm/controls/OrbitControls.js';

// HZ - BEGIN - 21/05/2024
/*La fonction movingAverage est utilisée pour lisser les données de vitesse dans le contexte de la 
 courbe de rotation galactique.
Explication de la Courbe de Rotation Galactique
Dans une simulation de galaxie, les étoiles tournent autour du centre galactique, 
et cette rotation est représentée graphiquement par une courbe de rotation galactique. 
Cette courbe montre comment la vitesse de rotation des étoiles varie en fonction de leur
distancce au centre de la galaxie.

1. Calcul de la vitesse totale
Au lieu de calculer uniquement la vitesse circulaire (qui est la vitesse le long d'une orbite circulaire), 
on calcule la norme absolue de la vitesse. Cela inclut toutes les composantes de la vitesse 
(radiale, tangente, etc.).
2. Division en portions
Tu divises le rayon de la galaxie en 50 segments. Pour chaque segment, on calcules la moyenne des 
vitesses des étoiles situées dans ce segment. Cela donne une idée de la distribution des vitesses 
à différentes distances du centre de la galaxie.
3. Moyenne mobile
Pour lisser les fluctuations des vitesses, on peu utiliser une moyenne mobile. Cela signifie que 
qu'on calcules la moyenne des vitesses sur plusieurs segments voisins. Par exemple, si on utilises 
une fenêtre de 5 segments, on prends la moyenne des vitesses des 5 segments autour de chaque point 
pour obtenir une courbe plus lisse.*/

// La fonction prend en parametre un tableau de nombres (les valeurs à lisser) et la taille de la 
// fenêtre pour la moyenne mobile (le nombre d'éléments à considérer pour chaque moyenne).
function movingAverage(data, windowSize) {
    // Initialiser un tableau vide pour stocker les résultats des moyennes mobiles
    let result = [];
    
    // Parcourir le tableau data avec une boucle
    // La boucle s'arrête pour s'assurer qu'il reste assez d'élément pour une fenêtre complète
    for (let i = 0; i < data.length - windowSize + 1; i++) {
        // Extraire une sous-section du tableau data de taille windowSize commencant à l'index i
        let window = data.slice(i, i + windowSize);
        
        // Calculer la somme des éléments dans la fenêtre actuelle
        // reduce utilisé pour accumuler la somme des éléments dans window
        let sum = window.reduce((acc, val) => acc + val, 0);
        
        // Calculer la moyenne de la fenêtre actuelle et l'ajouter au tableau result
        // La moyenne est obtenue en divisant la somme par windowSize
        result.push(sum / windowSize);
    }
    
    // Retourner le tableau result qui contient les moyennes mobiles calculées
    return result;
}
// HZ - END - 21/05/2024

class Cloud {
    constructor(stars) {
      this.stars = stars;
      this.geometry = null;
      this.points = null;
      this.createCloud();
    }
      createCloud() {
        const positions = new Float32Array(this.stars.length * 3);
        const colors = new Float32Array(stars.length * 3);
        this.stars.forEach((star, index) => {
          positions[index * 3] = star.posX;
          positions[index * 3 + 1] = star.posY;
          positions[index * 3 + 2] = star.posZ;
          if (star.massG > 0) {
            colors[index * 3 + 0] = 1;
            colors[index * 3 + 1] = 1;
            colors[index * 3 + 2] = 0;
          } else {
            colors[index * 3 + 0] = 1;
            colors[index * 3 + 1] = 0;
            colors[index * 3 + 2] = 1;
          }
        });
    
        const geometry = new THREE.BufferGeometry().setAttribute('position', new THREE.BufferAttribute(new Float32Array(positions), 3));
        geometry.setAttribute('color', new THREE.BufferAttribute(colors, 3));

        const material = new THREE.PointsMaterial({ size: 0.1,vertexColors: true , sizeAttenuation: true,blending: THREE.AdditiveBlending,transparent: true,opacity: 0.5 });
        this.geometry = geometry;
        this.points = new THREE.Points(geometry, material);
        scene.add(this.points);
      }
    
      updateCloud() {
        const positions = new Float32Array(this.stars.length * 3);
        this.stars.forEach((star, index) => {
          positions[index * 3] = star.posX;
          positions[index * 3 + 1] = star.posY;
          positions[index * 3 + 2] = star.posZ;
        });
      
        this.geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
      }
      
  
    removeCloud() {
      scene.remove(this.points);
      this.geometry.dispose();
      this.points.material.dispose();
      this.points = null;
    }
  }
  
  
class Star {
    constructor(massG, massI, posX, posY, posZ, velX, velY, velZ, color,type) {
        this.type = type;
        this.massG = massG; // Gravitational mass
        this.massI = massI; // Inertial mass
        this.posX = posX; // position x
        this.posY = posY; // position y
        this.posZ = posZ; // position y
        this.velX = velX; // vitesse x
        this.velY = velY; // vitesse y
        this.velZ = velZ; // vitesse y
        this.color = color; // Color of the star
        this.point = null; // point mesh 
    }
    AddToCloud(){

    }
    RemoveToCloud(){
        
    }
      createPoint() {
        let position = new THREE.Vector3(this.posX, this.posY, this.posZ);
        const geometry = new THREE.BufferGeometry().setAttribute('position', new THREE.BufferAttribute(new Float32Array([position.x, position.y, position.z]), 3));
        const material = new THREE.PointsMaterial({ color: this.color, size: 0.1, sizeAttenuation: true });
        this.point = new THREE.Points(geometry, material);
        this.point.position.set(0, 0, 0);
        scene.add(this.point);
      }
      
      removePoint(){

        scene.remove(this.point);
        this.point.geometry.dispose();
        this.point.material.dispose();
        this.point = null;
      }
    updatePointPositionReplay(starsReplayX,starsReplayY,starsReplayZ){
        this.point.geometry.attributes.position.setXYZ(0,starsReplayX,starsReplayY,starsReplayZ);
        this.point.geometry.attributes.position.needsUpdate = true;
    }
    updatePointPosition(){
        this.point.geometry.attributes.position.setXYZ(0, this.posX, this.posY, this.posZ);
        this.point.geometry.attributes.position.needsUpdate = true;
    }
    updatePositionReplay(starsReplayX,starsReplayY,starsReplayZ){ // starsReplay  position at timestep Replay of the particle
        this.posX = starsReplayX;
        this.posY = starsReplayY;
        this.posZ = starsReplayZ;
        this.updatePointPosition()
    }
      

    updatePosition(deltaTime, forceX, forceY, forceZ,BoundR) {

        const accelX = forceX / this.massI;
        const accelY = forceY / this.massI;
        const accelZ = forceZ / this.massI;


        this.posX += this.velX * deltaTime;
        this.posY += this.velY * deltaTime;
        this.posZ += this.velZ * deltaTime;


        this.velX += accelX * deltaTime;
        this.velY += accelY * deltaTime;
        this.velZ += accelZ * deltaTime;

        //this.checkBoundaryCondition(speed_init,angleVariance,BoundR);


        this.point.geometry.attributes.position.setXYZ(0, this.posX, this.posY, this.posZ);
        this.point.geometry.attributes.position.needsUpdate = true;

    }
    updatePositionLPReplay(deltaTime, forceX, forceY, forceZ,BoundR) {

        const accelX = forceX / this.massI;
        const accelY = forceY / this.massI;
        const accelZ = forceZ / this.massI;
        




        this.posX += this.velX * deltaTime + 1/2*accelX*deltaTime*deltaTime ;
        this.posY += this.velY * deltaTime + 1/2*accelY*deltaTime*deltaTime ;
        this.posZ += this.velZ * deltaTime + 1/2*accelZ*deltaTime*deltaTime ;

        this.checkBoundaryCondition(speed_init,angleVariance,BoundR);


        //this.point.geometry.attributes.position.setXYZ(0, this.posX, this.posY, this.posZ);
        //this.point.geometry.attributes.position.needsUpdate = true;

    }

    updatePositionLP(deltaTime, forceX, forceY, forceZ,BoundR) {

        const accelX = forceX / this.massI;
        const accelY = forceY / this.massI;
        const accelZ = forceZ / this.massI;
        




        this.posX += this.velX * deltaTime + 1/2*accelX*deltaTime*deltaTime ;
        this.posY += this.velY * deltaTime + 1/2*accelY*deltaTime*deltaTime ;
        this.posZ += this.velZ * deltaTime + 1/2*accelZ*deltaTime*deltaTime ;

        this.checkBoundaryCondition(speed_init,angleVariance,BoundR);

        // TEST ONE BUFFER...............

        //this.point.geometry.attributes.position.setXYZ(0, this.posX, this.posY, this.posZ);
        //this.point.geometry.attributes.position.needsUpdate = true;

    }
    updateVelocityLP(deltaTime,NewforceX, NewforceY, NewforceZ, forceX, forceY, forceZ) {

        const accelX = forceX / this.massI;
        const accelY = forceY / this.massI;
        const accelZ = forceZ / this.massI;

        const NewaccelX = NewforceX / this.massI;
        const NewaccelY = NewforceY / this.massI;
        const NewaccelZ = NewforceZ / this.massI;
        

        this.velX += 1/2*(accelX+NewaccelX) * deltaTime;
        this.velY += 1/2*(accelY+NewaccelY) * deltaTime;
        this.velZ += 1/2*(accelZ+NewaccelZ) * deltaTime;

    }
    updatePosition2(BoundR) {
        this.checkBoundaryCondition(speed_init,angleVariance,BoundR);
        this.point.geometry.attributes.position.setXYZ(0, this.posX, this.posY, this.posZ);
        this.point.geometry.attributes.position.needsUpdate = true;
       
    }
    checkBoundaryCondition(initialDarkMatterSpeed,angleVariance,BoundRadius) {
        const boundaryRadius = BoundRadius; // Half of the canvas size
        const distanceFromCenter = Math.sqrt(this.posX * this.posX  + this.posY*this.posY + this.posZ*this.posZ);
        

        if (distanceFromCenter > boundaryRadius) {
            //console.log('matter is out',distanceFromCenter,boundaryRadius)

            if (this.type === 'matter') {
                // Move the matter particle outside the canvas to cancel effect far away , -> useless condition , for auto gravitation matter never go outside
                /*this.posX =0
                this.posY =0
                this.posZ =boundaryRadius
                
                this.velX = 0;
                this.velY = 0;
                this.velZ = 0;*/

            } else if (this.type === 'darkMatter') {

                // Place the Nega matter particle at a random position along the circular boundary -> to have effect of "zero divergence" far the galaxy at boundary simulation
                // the size of galaxy need to be small to "dont affect"  Nega matter at boundary , and not have no desirable effet of flow NegaMatter


                let posX = 0
                let posY = 0
                let posZ =0

                

                let isgood2 = false;
                let distance=0
                while (isgood2 == false) {
                    posX = Math.random()  * boundaryRadius*2 - boundaryRadius;
                    posY = Math.random()  * boundaryRadius*2 - boundaryRadius;
                    posZ = Math.random()  * boundaryRadius*2 - boundaryRadius;
        
                     distance = Math.sqrt(posX * posX + posY * posY + posZ * posZ);
        
                    if (distance <= boundaryRadius) {
                        isgood2 = true
                    }
                }



                this.posX = ( posX/distance)*boundaryRadius
                this.posY = ( posY/distance)*boundaryRadius
                this.posZ = ( posZ/distance)*boundaryRadius



                   // Calculate the inward direction vector components
                let dx = - this.posX;
                let dy = - this.posY;
                let dz = - this.posZ;

                // Convert vector to spherical coordinates (r, theta, phi)
                let r = Math.sqrt(dx*dx + dy*dy + dz*dz);
                let theta2 = Math.acos(dz / r); // Polar angle from Z-axis
                let phi2 = Math.atan2(dy, dx); // Azimuthal angle from X-axis


                const maxThetaPerturbation = Math.PI / 4; // Adjust this value for less/more variance
                theta2 += (Math.random() - 0.5) * 2 * maxThetaPerturbation; // Randomly perturb theta within the allowed range


                theta2 = Math.max(0, Math.min(Math.PI, theta2));

                // Convert back to Cartesian coordinates
                dx = r * Math.sin(theta2) * Math.cos(phi2);
                dy = r * Math.sin(theta2) * Math.sin(phi2);
                dz = r * Math.cos(theta2);

                // Apply the speed to the direction vector components

                let isgood = false;
                let velX = 0;
                let velY = 0;
                let velZ = 0;
                let absVel=0

                while (isgood == false) {
                    velX = Math.random()  * initialDarkMatterSpeed*2 - initialDarkMatterSpeed;
                    velY = Math.random()  * initialDarkMatterSpeed*2 - initialDarkMatterSpeed;
                    velZ = Math.random()  * initialDarkMatterSpeed*2 - initialDarkMatterSpeed;

                    absVel = Math.sqrt(velX * velX + velY * velY + velZ * velZ);

                    if (absVel <= initialDarkMatterSpeed && absVel >= initialDarkMatterSpeed*0.8) {
                        isgood = true
                    }
                }
                this.velX = velX;
                this.velY = velY;
                this.velZ = velZ ;
            }
        }
    }

}


// Function to generate random values within a specified range
function randomBetween(min, max) {
    return Math.random() * (max - min) + min;
}

function initMatterCluster(numStars,densityNeg,Size,elipseX,elipseZ,speedscale,color='rgba(130, 255, 255, 0.3)') {
    const matterStars = [];
    for (let i = 0; i < numStars; i++) {

        const matterColor = color;

        const mass = densityNeg/numStars; // 10^12/10^4 ( total galaxy solar mass / nmbr star)
        const angle = randomBetween(0, 2 * Math.PI);
        const radius = randomBetween(0.01, Size); // Distance from the center (in pixel)
        const posX =  radius * Math.cos(angle) * elipseX;
        const posY =  radius * Math.sin(angle) ;
        const posZ = randomBetween(0.01, Size/4);


        let ad=speedscale;

        // Velocity for rotational motion
      
        let velX = 0;
        let velY = 0;
        let velZ = 0;

        if ( radius < Size/5 ){
            velX = -Math.sin(angle) *( ad*radius+200);
            velY = Math.cos(angle) * (ad*radius+200);
            velZ = 0
        }
        else{
             velX = -Math.sin(angle) * ad*radius;
             velY = Math.cos(angle) * ad*radius;
             velZ = 0
        }
        
        

        matterStars.push(new Star(mass, Math.abs(mass), posX, posY, posZ, velX, velY,velZ,matterColor,'matter'));
    }
    return matterStars;
}


function initDarkMatterCluster(numStars, speed,densityNeg,size,simulationRadius,color='rgba(255, 0, 0, 1)') {
    const darkMatterStars = [];
    const darkMatterColor =color; 


    for (let i = 0; i < numStars; i++) {
        //const mass = densityNeg/numStars;
        const mass =densityNeg*4/3*Math.PI*simulationRadius*simulationRadius*simulationRadius/numStars

      
        let posX=0 ;
        let posY=0;
        let posZ =0;
        let distance=0;
        let isgood = false;

        while (isgood == false) {
             posX = Math.random()  * simulationRadius*2 - simulationRadius;
             posY = Math.random()  * simulationRadius*2 - simulationRadius;
             posZ = Math.random()  * simulationRadius*2 - simulationRadius;

             distance = Math.sqrt(posX * posX + posY * posY + posZ * posZ);

            if (distance <= simulationRadius && distance >= size) {
                isgood = true
            }
        }


        distance = Math.sqrt(posX * posX + posY * posY + posZ * posZ);

        isgood = false;
        let velX = 0;
        let velY = 0;
        let velZ = 0;
        let absVel=0


     //try to mimic pseudo uniforme random direc
        while (isgood == false) {
            velX = Math.random()  * speed*2 - speed;
            velY = Math.random()  * speed*2 - speed;
            velZ = Math.random()  * speed*2 - speed;

            absVel = Math.sqrt(velX * velX + velY * velY + velZ * velZ);

            if (absVel <= speed && absVel >= speed*0.8) {
                isgood = true
            }
        }
        //console.log('dark position init',radi)
        

        
       // console.log(posX,posY,posZ,velX,velY,velZ)
        darkMatterStars.push(new Star(mass, Math.abs(mass), posX, posY,posZ, velX, velY,velZ,darkMatterColor,'darkMatter'));
    }
    return darkMatterStars;
}

let computeForcesKernel;
let computeForcesKernel2;

function initializeComputeForcesKernel(numStars) {
    // If a kernel already exists, destroy it
    if (computeForcesKernel) {
        computeForcesKernel.destroy();
    }
    if (computeForcesKernel2) {
        computeForcesKernel2.destroy();
    }
    // Initialize GPU instance if not already initialized or if destroyed
    console.log('init gpu')
    // HZ - BEGIN - 21/05/2024
    //const gpu = new GPU.GPU({mode: 'gpu'});
    const gpu = new GPU();
    // HZ - END - 21/05/2024
    gpu.maxLoopSize = 30000000
    // Adjust the output size based on the number of stars
    const loopSize = numStars/2; //starposi + starnega
    const outputSizeChunked = numStars/2 ;
    // Create a new kernel with the updated output size
    console.log('init kernel')
    computeForcesKernel = gpu.createKernel(function(starMasses, starPositionsX, starPositionsY,starPositionsZ,starMasses2, starPositionsX2, starPositionsY2,starPositionsZ2, G, densityNeg) {
        let forceX = 0;
        let forceY = 0;
        let forceZ = 0;
 
    
        const myPosX = starPositionsX[this.thread.x];
        const myPosY = starPositionsY[this.thread.x];
        const myPosZ = starPositionsZ[this.thread.x];
        const myMass = starMasses[this.thread.x];

    
       for (let i = 0 ; i < this.constants.loop ; i++) {
            if (i !== (this.thread.x )) {
                const dx = starPositionsX2[i] - myPosX;
                const dy = starPositionsY2[i] - myPosY;
                const dz = starPositionsZ2[i] - myPosZ;
                const cutoff = 2; //smouthing gravity short scale epsilon
                let distance= Math.sqrt(dx * dx + dy * dy + dz * dz );

                if(distance >= cutoff) {
                    const forceMagnitude = G * myMass * starMasses2[i] / (distance * distance);
                    //f1ratio += forceMagnitude
                    forceX += forceMagnitude * (dx / distance);
                    forceY += forceMagnitude * (dy / distance);
                    forceZ += forceMagnitude * (dz / distance);
                }

            }
        }
    
        // Additional fictive force due to constant Negative matter density outside (homogene density)

        const distanceFromCentere = Math.sqrt(myPosX  * myPosX  + myPosY * myPosY + myPosZ * myPosZ);
        //console.log("distanceFromCenter", distanceFromCentere)
        

        // fictive "shadow of negaMatter outside the simulation" (Big lacune + Nega matter inside = 0 (infinite homogeneous repartition) ) 
        const forceMagnitude = 4/3*Math.PI*G*myMass*densityNeg*distanceFromCentere / 2 // interior of the constante density Shadow NegaBigLacune of simulation, divide by 2 because of the two chunk that pass 2 time on this force



        forceX += forceMagnitude * ((myPosX ) / distanceFromCentere);
        forceY += forceMagnitude * ((myPosY ) / distanceFromCentere);
        forceZ += forceMagnitude * ((myPosZ ) / distanceFromCentere);

        

        //const ratio = f1ratio/f2ratio;
        return [forceX, forceY,forceZ];
    }, {
        constants: { loop : loopSize },
        output: [ outputSizeChunked ],
    });
    computeForcesKernel.setLoopMaxIterations(gpu.maxLoopSize);



    computeForcesKernel2 = gpu.createKernel(function(force1,force2) {

        let forceTotX = force1[this.thread.x][0] + force2[this.thread.x][0]

        let forceTotY = force1[this.thread.x][1] + force2[this.thread.x][1]

        let forceTotZ = force1[this.thread.x][2] + force2[this.thread.x][2]

        return [forceTotX,forceTotY,forceTotZ];
    }, {

        output: [ outputSizeChunked ],
    });

}





let isFisrtLoop=true

let forceNow =[]
let forceNew=[]
let forceNowX=[]
let forceNowY=[]
let forceNowZ=[]

function updateGalaxyGPULP(stars, deltaTime, G, densityNeg, boundaryRadius) {
    //console.log("updateGPU")
    // Prepare data for GPU
    
    for (let i = 0; i < stars.length; i++) {

        stars[i].updatePositionLP(deltaTime, forceNowX[i], forceNowY[i],forceNowZ[i],boundaryRadius);
 
     }
     //console.log("star posi",stars)



    const starMasses = stars.map(star => star.massG);
    const starPositionsX = stars.map(star => star.posX);
    const starPositionsY = stars.map(star => star.posY);
    const starPositionsZ = stars.map(star => star.posZ);

    // First kernel to compute forces
    let chunckN=0 ///////////////
    let forcesX2 = []
    let forcesY2 = []
    let forcesZ2 = []
  


    

        const halfLength = Math.floor(starMasses.length / 2);

        let SubMasses = starMasses.slice(0, halfLength);
        let SubPosX = starPositionsX.slice(0, halfLength);
        let SubPosY = starPositionsY.slice(0, halfLength);
        let SubPosZ = starPositionsZ.slice(0, halfLength);

        let SubMasses2 = starMasses.slice(halfLength,-1);
        let SubPosX2 = starPositionsX.slice(halfLength,-1);
        let SubPosY2 = starPositionsY.slice(halfLength,-1);
        let SubPosZ2 = starPositionsZ.slice(halfLength,-1);

        let forcesT = computeForcesKernel(SubMasses, SubPosX, SubPosY, SubPosZ,SubMasses, SubPosX, SubPosY, SubPosZ, G, densityNeg);
        //console.log("forcesT", forcesT);


        let forcesT2 = computeForcesKernel(SubMasses, SubPosX, SubPosY, SubPosZ,SubMasses2, SubPosX2, SubPosY2, SubPosZ2, G, densityNeg);

        //console.log("forcesT2", forcesT2);

        let forcetot = computeForcesKernel2(forcesT,forcesT2)

        //console.log("forcetot", forcetot);

       
    
        let forXx = forcetot.map(force => force[0]);
        let foryy = forcetot.map(force => force[1]);
        let forzz = forcetot.map(force => force[2]);
    
        forcesX2 = forcesX2.concat(forXx);
        forcesY2 =  forcesY2.concat(foryy);
        forcesZ2 = forcesZ2.concat(forzz);

    
         forcesT = computeForcesKernel(SubMasses2, SubPosX2, SubPosY2, SubPosZ2,SubMasses, SubPosX, SubPosY, SubPosZ, G, densityNeg);

         forcesT2 = computeForcesKernel(SubMasses2, SubPosX2, SubPosY2, SubPosZ2,SubMasses2, SubPosX2, SubPosY2, SubPosZ2, G, densityNeg);
        
         forcetot = computeForcesKernel2(forcesT,forcesT2)

         let forXx2 = forcetot.map(force => force[0]);
         let foryy2 = forcetot.map(force => force[1]);
         let forzz2 = forcetot.map(force => force[2]);
    
         forcesX2 = forcesX2.concat(forXx2);
         forcesY2 =  forcesY2.concat(foryy2);
         forcesZ2 = forcesZ2.concat(forzz2);




        
        
        
    

   // console.log("forcesX", forcesX2);
    //console.log("forcesY", forcesY2);
    //console.log("forcesZ", forcesZ2);

    
    // Update the stars array with new positions and velocities cpu

    for (let i = 0; i < stars.length; i++) {

       stars[i].updateVelocityLP(deltaTime, forcesX2[i], forcesY2[i],forcesZ2[i],forceNowX[i], forceNowY[i],forceNowZ[i],boundaryRadius);

    }
    
    forceNowX = forcesX2;
    forceNowY = forcesY2 
    forceNowZ = forcesZ2
    
}


let ReplayArray=[]
let VelocityReplay=[]

const MAX_ENTRIES = 10;  // Adjust based on your empirical memory usage testing
let entries_index=0
let fileIndex = 0;
let directoryHandle = null;
const isFileSystemAccessSupported = 'showSaveFilePicker' in window;

async function updateGalaxyGPULPReplay(stars, deltaTime, G, densityNeg, boundaryRadius) {
    // Prepare data for GPU
    let stepReplayArray=[]
    for (let i = 0; i < stars.length; i++) {

        stars[i].updatePositionLPReplay(deltaTime, forceNowX[i], forceNowY[i],forceNowZ[i],boundaryRadius);

        let arrayStarPos = [stars[i].posX,stars[i].posY,stars[i].posZ]
        stepReplayArray.push(arrayStarPos)
     }
     ReplayArray.push(stepReplayArray)

     //console.log("ReplayArray",ReplayArray)



    const starMasses = stars.map(star => star.massG);
    const starPositionsX = stars.map(star => star.posX);
    const starPositionsY = stars.map(star => star.posY);
    const starPositionsZ = stars.map(star => star.posZ);

    // First kernel to compute forces

    /*const forces = computeForcesKernel(starMasses, starPositionsX, starPositionsY,starPositionsZ, G, densityNeg,2);
    forceNew = forces
    const forcesX = forces.map(force => force[0]);
    const forcesY = forces.map(force => force[1]);
    const forcesZ = forces.map(force => force[2]);*/

    
    // First kernel to compute forces
    /*let chunckN=0 ///////////////
    let forcesX = []
    let forcesY = []
    let forcesZ = []


    for (let i = 0; i < 2; i++) {
        let forces = computeForcesKernel(starMasses, starPositionsX, starPositionsY, starPositionsZ, G, densityNeg, chunckN);
        
        chunckN++;
    
        let forXx = forces.map(force => force[0]);
        let foryy = forces.map(force => force[1]);
        let forzz = forces.map(force => force[2]);
    
        forcesX= forcesX.concat(forXx);
        forcesY =  forcesY.concat(foryy);
        forcesZ = forcesZ.concat(forzz);
    }
    //console.log("forcesX", forcesX);*/
    let forcesX = []
    let forcesY = []
    let forcesZ = []
    
    const halfLength = Math.floor(starMasses.length / 2);

    let SubMasses = starMasses.slice(0, halfLength);
    let SubPosX = starPositionsX.slice(0, halfLength);
    let SubPosY = starPositionsY.slice(0, halfLength);
    let SubPosZ = starPositionsZ.slice(0, halfLength);

    let SubMasses2 = starMasses.slice(halfLength);
    let SubPosX2 = starPositionsX.slice(halfLength);
    let SubPosY2 = starPositionsY.slice(halfLength);
    let SubPosZ2 = starPositionsZ.slice(halfLength);

    let forcesT = computeForcesKernel(SubMasses, SubPosX, SubPosY, SubPosZ,SubMasses, SubPosX, SubPosY, SubPosZ, G, densityNeg);
    //console.log("forcesT", forcesT);


    let forcesT2 = computeForcesKernel(SubMasses, SubPosX, SubPosY, SubPosZ,SubMasses2, SubPosX2, SubPosY2, SubPosZ2, G, densityNeg);

    //console.log("forcesT2", forcesT2);

    let forcetot = computeForcesKernel2(forcesT,forcesT2)

    //console.log("forcetot", forcetot);

   

    let forXx = forcetot.map(force => force[0]);
    let foryy = forcetot.map(force => force[1]);
    let forzz = forcetot.map(force => force[2]);

    forcesX = forcesX.concat(forXx);
    forcesY =  forcesY.concat(foryy);
    forcesZ = forcesZ.concat(forzz);


     forcesT = computeForcesKernel(SubMasses2, SubPosX2, SubPosY2, SubPosZ2,SubMasses, SubPosX, SubPosY, SubPosZ, G, densityNeg);

     forcesT2 = computeForcesKernel(SubMasses2, SubPosX2, SubPosY2, SubPosZ2,SubMasses2, SubPosX2, SubPosY2, SubPosZ2, G, densityNeg);
    
     forcetot = computeForcesKernel2(forcesT,forcesT2)


     let forXx2 = forcetot.map(force => force[0]);
     let foryy2 = forcetot.map(force => force[1]);
     let forzz2 = forcetot.map(force => force[2]);

     forcesX = forcesX.concat(forXx2);
     forcesY =  forcesY.concat(foryy2);
     forcesZ = forcesZ.concat(forzz2);

    


    
    // Update the stars array with new positions and velocities cpu
    let stepReplayVelocity=[]

    for (let i = 0; i < stars.length; i++) {

       stars[i].updateVelocityLP(deltaTime, forcesX[i], forcesY[i],forcesZ[i],forceNowX[i], forceNowY[i],forceNowZ[i],boundaryRadius);
       let arrayStarVel = [stars[i].velX,stars[i].velY,stars[i].velZ]
       stepReplayVelocity.push(arrayStarVel)
    }
    VelocityReplay.push(stepReplayVelocity)
    //console.log("VelocityReplay",VelocityReplay)

    forceNowX = forcesX;
    forceNowY = forcesY;
    forceNowZ = forcesZ;
    entries_index++

    // Check if we need to save to file and clear arrays
    if (stepNow % MAX_ENTRIES == 0 && stepNow!=0 ) {
        //console.log("enter in the condition saving for max entry 10")
        entries_index=0
        saveDataAndClearArrays(function() {
            // After save callback
            //console.log('Data saved and arrays cleared.');
        });
         //ReplayArray=[]
         //VelocityReplay=[]

        
    }
}
/*document.getElementById('selectFolderButton').addEventListener('click', async function() {
    if (isFileSystemAccessSupported) {
        await selectFolder();
        console.log('Folder selected.');
    } else {
        alert('Your browser does not support the File System Access API');
    }
});*/

async function selectFolder() {
    try {
        directoryHandle = await window.showDirectoryPicker();
        console.log('Folder selected:', directoryHandle);
    } catch (error) {
        console.error('Directory selection failed:', error);
    }
}

async function saveDataAndClearArrays() {
    if (ReplayArray.length !== VelocityReplay.length) {
        console.error('Data inconsistency detected before saving. ReplayArray length:', ReplayArray.length, 'VelocityReplay length:', VelocityReplay.length);
        return;
    }

    const csvDataPos = convertToCSV(ReplayArray.flat());
    const csvDataVel = convertToCSV(VelocityReplay.flat());

    if (directoryHandle) {
        try {
            await appendToFile(directoryHandle, csvDataPos, `replay_positions.csv`);
            await appendToFile(directoryHandle, csvDataVel, `replay_velocities.csv`);
            // Clear arrays after successful file operations
            ReplayArray = [];
            VelocityReplay = [];
        } catch (error) {
            console.error('Error during file write operation:', error);
            throw error; // Rethrow the error to handle it in the calling function
        }
    } else {
        console.warn('No directory selected. Saving to default download folder.');
        appendToFileDefault(csvDataPos, `replay_positions.csv`);
        appendToFileDefault(csvDataVel, `replay_velocities.csv`);
        // Clear arrays after default file operations
        ReplayArray = [];
        VelocityReplay = [];
    }

    await validateSavedData();
}

async function appendToFile(directoryHandle, data, filename) {
    try {
        const fileHandle = await directoryHandle.getFileHandle(filename, { create: true });
        const writable = await fileHandle.createWritable({ keepExistingData: true });
        await writable.write(data);
        await writable.close();
    } catch (error) {
        console.error('Error during file write operation:', error);
        throw error; // Rethrow the error to be caught in saveDataAndClearArrays
    }
}

function appendToFileDefault(data, filename) {
    const existingContent = localStorage.getItem(filename) || '';
    const newContent = existingContent + data + '\n';
    localStorage.setItem(filename, newContent);

    const blob = new Blob([newContent], { type: 'text/csv' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
}

function convertToCSV(dataArray) {
    return dataArray.map(row => row.join(',')).join('\n');
}

async function validateSavedData() {
    console.log("validateSavedData enter function")
    if (directoryHandle) {
        const dataPos = await readFile(directoryHandle, 'replay_positions.csv');
        const dataVel = await readFile(directoryHandle, 'replay_velocities.csv');
        const posLines = dataPos.split('\n').length;
        const velLines = dataVel.split('\n').length;
        console.log('Validation: Position data lines:', posLines);
        console.log('Validation: Velocity data lines:', velLines);
        if (posLines !== velLines) {
            console.error('Data inconsistency detected after saving. Position data lines:', posLines, 'Velocity data lines:', velLines);
        } else {
            console.log('Data consistency verified after saving.');
        }
    } else {
        console.log('Validation skipped, no directory handle available.');
    }
}

async function readFile(directoryHandle, filename) {
    try {
        const fileHandle = await directoryHandle.getFileHandle(filename);
        const file = await fileHandle.getFile();
        return await file.text();
    } catch (error) {
        console.error('File read failed:', error);
        throw error;
    }
}
// Create a scene, camera, and renderer
let scene = new THREE.Scene();
let camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 10000);
let renderer = new THREE.WebGLRenderer();
renderer.setSize(window.innerWidth, window.innerHeight);
document.body.appendChild(renderer.domElement);

// Create an OrbitControls object
let controls = new OrbitControls(camera, renderer.domElement);

// Set the position and target of the camera
camera.position.set(0, 0, 800);
controls.target.set(0, 0, 0);
controls.update();

let Radius_SImulation= 1600;
let scale= 0.0005;
let  deltaTime= 1000000;
let  numStarsNeg= 1000;
let   numStarsPos=1000;
let densityNeg= -100000000000;
let TotmassPos= 100000000000;
let hole= 20;
let galactR= 90;
let elips= 1;
let speedRot= 0.00000001;
let speed_init= 0.0000000000000000001;
let colorNeg= '#ff0000';
let colorPos= '#00ff00';
let G= 0;
let timestep = 200;

const angleVariance = Math.PI ;

let stars = [];
let cloud = null;
let animationFrameId;
let animationFrameId2;
let animationFrameId3;
let animationFrameId4;

let darkMatterDensity=0;

var params = {
    Radius_SImulation: 500,
    scale: 0.0005,
    deltaTime: 0.002,
    numStarsNeg: 100000,
    numStarsPos:20000,
    densityNeg: -0.1*Math.pow(10,5), // solar mass / ly^3 (density mass of nega stars)
    TotmassPos: 5*Math.pow(10,11), // ~ milky way solar mass (total mass of the galaxy)
    hole: 35,
    galactR: 25, // 25 unit x 1/scale = 50 000 ly ~ milky way radius
    elips: 1,
    speedRot: 30, //  rad per billion years (10^9 year)
    speed_init: 500, // ly/Gy
    colorNeg: '#ff0000',
    colorPos: '#00ff00',
    timestep:100,
    initialize: function() {
        Radius_SImulation= params.Radius_SImulation;
        scale = params.scale;
        deltaTime= params.deltaTime;
        numStarsNeg= params.numStarsNeg;
        numStarsPos=params.numStarsPos;
        densityNeg= params.densityNeg;
        TotmassPos= params.TotmassPos;
        hole= params.hole;
        galactR= params.galactR;
        elips= params.elips;
        speedRot= params.speedRot;
        speed_init= params.speed_init;
        colorNeg=params.colorNeg;
        colorPos=params.colorPos;
        G = params.G;
        timestep =params.timestep;

        // Code to initialize the simulation
        //G = 1.56*Math.pow(10,-13) * scale * scale * scale  // Adjusted G from ly^3 / (M☉ * yr^2) to px^3 / (M☉ * yr^2) 

        G = 155800 * scale * scale * scale // unit_engine^3 / (M☉ * Gyr^2) G = 155800 -> ly^3 / (M☉ * Gyr^2) Gyr=10^9 yr
        //console.log("params.G",G)
        let totalstars = numStarsNeg +  numStarsPos
        console.log("G ",G )
        initializeComputeForcesKernel(totalstars);
        initSimulation()
    },
    runRealTime: function() {
        // Code to run the simulation in real time
        //animate2();
        //animateLeap()
        animateLp();
    },
    runWait: function() {
        RunSimuWorker(timestep)
    },
    replay: function() {
        Replay(timestep)
    }
};

var gui = new dat.GUI();


gui.add(params, 'Radius_SImulation', hole+1, 10000).step(10).name('Radius of Simulation in light year').onChange(function(value) {
    Radius_SImulation = parseFloat(value,10);
});

gui.add(params, 'scale', 0, 1).step(0.0001).name('scale 1 unit / 2000 ly = 0.0005').onChange(function(value) {
    scale = parseFloat(value,10);
    
});

gui.add(params, 'deltaTime', 0, 1).step(0.00001).name('deltaTime in Giga-year').onChange(function(value) {
    deltaTime = parseInt(value, 10);
});
gui.add(params, 'timestep', 0, 100000).step(1).name('timestep').onChange(function(value) {
    timestep = parseFloat(value,10);
    
});

gui.add(params, 'numStarsNeg', 0, 10000000).step(1).name('number nega stars').onChange(function(value) {
    numStarsNeg = parseInt(value, 10);
});
gui.add(params, 'numStarsPos', 0, 10000000).step(1).name('number posi stars').onChange(function(value) {
    numStarsPos = parseInt(value, 10);
});

gui.add(params, 'densityNeg', -100000000000, 100000000000).step(0.1).name('density massNega solar mass / ly^3').onChange(function(value) {
    densityNeg = parseFloat(value,10);
});

gui.add(params, 'TotmassPos',  -10000000000000000, 10000000000000000).step(1000).name('total solar mass Galaxy').onChange(function(value) {
    TotmassPos = parseFloat(value,10);
});

gui.add(params, 'hole', 0, 1000).step(1).name('holei').onChange(function(value) {
    hole = parseInt(value, 10);
});

gui.add(params, 'galactR', 0, Radius_SImulation-1).step(1).name('galactR').onChange(function(value) {
    galactR = parseInt(value, 10);
});

gui.add(params, 'elips', 0, 1).step(0.1).name('elips').onChange(function(value) {
    elips = parseFloat(value,10);
});

gui.add(params, 'speedRot', 0, 1000).step(0.1).name('speedRot rad/Gy').onChange(function(value) {
    speedRot = parseFloat(value,10);
    
    
});

gui.add(params, 'speed_init', 0, 100000).step(1).name('speed_init ly/Gy').onChange(function(value) {
    speed_init = parseFloat(value,10);
    
});

// For colors, you'll need to use a color picker
/*gui.addColor(params, 'colorNeg').name('colorNeg').onChange(function(value) {
    colorNeg = value;
});*/

/*gui.addColor(params, 'colorPos').name('colorPos').onChange(function(value) {
    colorPos = value;
});*/

gui.add(params, 'initialize');
gui.add(params, 'runRealTime');
//gui.add(params, 'runWait');
//gui.add(params, 'replay');


let rotationCurveLine = null;

function plotGalacticRotationCurve(stars, numSegments = 50) {
    // Filtrer les étoiles de type matière ordinaire
    const matterStars = stars.filter(star => star.type === 'matter');

    // Initialiser les variables pour la masse totale et le centre de masse
    let totalMass = 0;
    let comX = 0, comY = 0, comZ = 0;

    // Calculer la masse totale et le centre de masse en sommant les positions pondérées par la masse
    matterStars.forEach(star => {
        totalMass += star.massG;
        comX += star.posX * star.massG;
        comY += star.posY * star.massG;
        comZ += star.posZ * star.massG;
    });

    // Calculer les coordonnées du centre de masse
    comX /= totalMass;
    comY /= totalMass;
    comZ /= totalMass;

    // Créer un tableau des rayons et des vitesses pour chaque étoile
    const data = matterStars.map(star => {
        const dx = star.posX - comX;
        const dz = star.posZ - comZ;
        const radius = Math.sqrt(dx * dx + dz * dz); // Rayon dans le plan x-z
        const velocityNorm = Math.sqrt(star.velX * star.velX + star.velY * star.velY + star.velZ * star.velZ); // Norme de la vitesse
        return { radius, velocityNorm };
    });

    // Trouver la vitesse maximale pour normaliser les donnée
    const maxVelocity = Math.max(...data.map(point => point.velocityNorm));

    // Trier les données par rayon croissant
    data.sort((a, b) => a.radius - b.radius);

    // Déterminer le rayon maximal et la taille des segments
    const maxRadius = data[data.length - 1].radius;
    const segmentSize = maxRadius / numSegments;

    // Initialiser les segments de donnée
    const segmentData = Array.from({ length: numSegments }, (_, i) => ({ radius: 0, avgVelocity: 0, count: 0 }));

    // Remplir les segments avec les données des étoiles
    data.forEach(point => {
        const segmentIndex = Math.floor(point.radius / segmentSize);
        if (segmentIndex < numSegments) {
            segmentData[segmentIndex].radius += point.radius;
            segmentData[segmentIndex].avgVelocity += point.velocityNorm;
            segmentData[segmentIndex].count += 1;
        }
    });

    // Calculer les valeurs moyennes pour chaque segment
    segmentData.forEach(segment => {
        if (segment.count > 0) {
            segment.radius /= segment.count;
            segment.avgVelocity /= segment.count;
        }
    });

    // Filtrer les segments vides
    const filteredSegmentData = segmentData.filter(segment => segment.count > 0);

    // HZ - BEGIN - 21/05/2024
    // Appliquer une moyenne mobile aux vitesses moyennes
    const avgVelocities = filteredSegmentData.map(segment => segment.avgVelocity);
    const smoothedVelocities = movingAverage(avgVelocities, 5); // Moyenne mobile avec une fenêtre de 5
    // HZ - END - 21/05/2024
    // Effacer le canvas pour dessiner la nouvelle courbe
    overlayContext.clearRect(0, 0, overlayCanvas.width, overlayCanvas.height);
    overlayContext.beginPath();
    overlayContext.moveTo(0, overlayCanvas.height);

    // HZ - BEGIN - 21/05/2024
    // Dessiner la courbe lissée
    smoothedVelocities.forEach((velocity, index) => {
        const x = (filteredSegmentData[index].radius / maxRadius) * overlayCanvas.width;
        const y = overlayCanvas.height - (velocity / maxVelocity * overlayCanvas.height);
        overlayContext.lineTo(x, y);
    });
    // HZ - END - 21/05/2024

    // Appliquer les styles et tracer la courbe
    overlayContext.strokeStyle = 'red';
    overlayContext.lineWidth = 2;
    overlayContext.stroke();

    // HZ - BEGIN - 21/05/2024
    // Ajouter les libellés des axes
    overlayContext.fillStyle = 'white';
    overlayContext.font = '16px Arial';
    // Libellé de l'axe des abscisses 
    overlayContext.fillText('Distance (kpc)', overlayCanvas.width / 2, overlayCanvas.height - 10);
    // Libellé de l'axe des ordonnées
    overlayContext.save();
    overlayContext.rotate(-Math.PI / 2);
    overlayContext.fillText('Vitesse Circulaire (km/s)', -overlayCanvas.height / 2, 20);
    overlayContext.restore();
    // HZ - END - 21/05/2024
}

let overlayCanvas, overlayContext;

// Initialize the overlay canvas
overlayCanvas = document.getElementById('overlayCanvas');
overlayContext = overlayCanvas.getContext('2d');

function initSimulation() {
    ReplayArray = []
    VelocityReplay=[];
    stepNow=0
    //initializeComputeForcesKernel(params.numStarsPos+params.numStarsNeg);

    

    // If there's an existing animation frame, cancel it
    if (animationFrameId) {
        cancelAnimationFrame(animationFrameId);
    }
    if (animationFrameId2) {
        cancelAnimationFrame(animationFrameId2);
    }
    if (animationFrameId3) {
        cancelAnimationFrame(animationFrameId3);
    }
    if (animationFrameId4) {
        cancelAnimationFrame(animationFrameId4);
    }
    
// TEST ONE BUFFER...................
    /*if(stars.length !=0){
        stars.forEach(star => {
        star.removePoint()
        
    })}*/

    if (cloud != null) {
        cloud.removeCloud();
      }
    scene.traverse((object) => {
    if (object.isMesh) {
        object.geometry.dispose();
        object.material.dispose();
    }
    });
    //renderer.dispose();
   //renderer.forceContextLoss();

    stars = [];
        // Draw the circular boundary
    //drawBoundary();
    // Assume initMatterCluster and initDarkMatterCluster are defined elsewhere
    const matterStars = initMatterCluster(numStarsPos, TotmassPos, galactR, elips,elips, speedRot,colorPos);
    //let speed_init = Math.sqrt(kB * temprature_init);
    const darkMatterStars = initDarkMatterCluster(numStarsNeg, speed_init, densityNeg, hole,Radius_SImulation,colorNeg );
    //darkMatterDensity = calculateDarkMatterDensity(darkMatterStars, Radius_SImulation);

    //densityNeg
    // Combine the stars
    stars = [...matterStars, ...darkMatterStars];
    cloud = new Cloud(stars);
    //cloud.setColorNeg(colorNeg);
    //cloud.setColorPos(colorPos)
    //cloud.setColorNeg(colorNeg);
    //console.log("stars.length",stars.length)

    // You can also draw the initial state of stars here if needed


    // TEST ONE BUFFER................................................................

    /*stars.forEach(star => {
        //star.draw();
        star.createPoint()

        //console.log("star.point",star.point)
    
    });*/

    const starMasses = stars.map(star => star.massG);
    const starPositionsX = stars.map(star => star.posX);
    const starPositionsY = stars.map(star => star.posY);
    const starPositionsZ = stars.map(star => star.posZ);
    console.log("first compute kernel calcul on init : START")
    //console.log("stars",stars)
    //forceNow = computeForcesKernel(starMasses, starPositionsX, starPositionsY,starPositionsZ, G, densityNeg)


    let chunckN2=0
    let forcesX22 = []
    let forcesY22 = []
    let forcesZ22 = []












        /*console.log(chunckN2)
        let forceNow1 = computeForcesKernel(starMasses, starPositionsX, starPositionsY,starPositionsZ, G, densityNeg,chunckN2,0)
        let forceNow2 = computeForcesKernel(starMasses, starPositionsX, starPositionsY,starPositionsZ, G, densityNeg,chunckN2,1)
        chunckN2++
        forceNow = computeForcesKernel2(forceNow1,forceNow2)

       // console.log("forceNow",forceNow)
        forceNew = forceNow

        let forX= forceNow.map(forceNow => forceNow[0])
        let fory= forceNow.map(forceNow => forceNow[1])
        let forz= forceNow.map(forceNow => forceNow[2])
        
        forcesX22=forcesX22.concat(forX);
        forcesY22=forcesY22.concat(fory);
        forcesZ22=forcesZ22.concat(forz);*/



        ///////////////////////////////// test Chunk the GPU calcul -> pas top ,permet +100K avant lost context gpu.js ( 2D array au lieu de 1D array ?)

        const halfLength = Math.floor(starMasses.length / 2);

        let SubMasses = starMasses.slice(0, halfLength);
        let SubPosX = starPositionsX.slice(0, halfLength);
        let SubPosY = starPositionsY.slice(0, halfLength);
        let SubPosZ = starPositionsZ.slice(0, halfLength);

        let SubMasses2 = starMasses.slice(halfLength);
        let SubPosX2 = starPositionsX.slice(halfLength);
        let SubPosY2 = starPositionsY.slice(halfLength);
        let SubPosZ2 = starPositionsZ.slice(halfLength);
        
        //console.log("iciii",starMasses.length )

        let forcesT = computeForcesKernel(SubMasses, SubPosX, SubPosY, SubPosZ,SubMasses, SubPosX, SubPosY, SubPosZ, G, densityNeg);



        let forcesT2 = computeForcesKernel(SubMasses, SubPosX, SubPosY, SubPosZ,SubMasses2, SubPosX2, SubPosY2, SubPosZ2, G, densityNeg);

        
        let forcetot= computeForcesKernel2(forcesT,forcesT2)


       
    
        let forXx = forcetot.map(force => force[0]);
        let foryy = forcetot.map(force => force[1]);
        let forzz = forcetot.map(force => force[2]);
    
        forcesX22 = forcesX22.concat(forXx);
        forcesY22 =  forcesY22.concat(foryy);
        forcesZ22 = forcesZ22.concat(forzz);

    
         forcesT = computeForcesKernel(SubMasses2, SubPosX2, SubPosY2, SubPosZ2, SubMasses, SubPosX, SubPosY, SubPosZ, G, densityNeg);

         forcesT2 = computeForcesKernel(SubMasses2, SubPosX2, SubPosY2, SubPosZ2,SubMasses2, SubPosX2, SubPosY2, SubPosZ2, G, densityNeg);
        
         forcetot= computeForcesKernel2(forcesT,forcesT2)

         let forXx2 = forcetot.map(force => force[0]);
         let foryy2 = forcetot.map(force => force[1]);
         let forzz2 = forcetot.map(force => force[2]);
    
         forcesX22 = forcesX22.concat(forXx2);
         forcesY22 =  forcesY22.concat(foryy2);
         forcesZ22 = forcesZ22.concat(forzz2);


 
     

    forceNowX = forcesX22;
    forceNowY = forcesY22;
    forceNowZ = forcesZ22;

    //console.log(forceNow)
    //console.log(forceNow)
    /*forceNowX = forceNow.map(forceNow => forceNow[0]);
    forceNowY = forceNow.map(forceNow => forceNow[1]);
    forceNowZ = forceNow.map(forceNow => forceNow[2]);*/
    //console.log("first RENDERER: START")
    renderer.render(scene, camera);

    plotGalacticRotationCurve(stars);
    //console.log("first RENDERER: END")
    animateFree()
}

let forcesTot=[]


function animateLp() {
    if (animationFrameId) {
        cancelAnimationFrame(animationFrameId);
    }
    if (animationFrameId2) {
        cancelAnimationFrame(animationFrameId2);
    }
    if (animationFrameId4) {
        cancelAnimationFrame(animationFrameId4);
    }
   
    animationFrameId3=requestAnimationFrame(animateLp);


    updateGalaxyGPULP(stars, deltaTime, G, densityNeg, Radius_SImulation);
    cloud.updateCloud();
    controls.update();
    plotGalacticRotationCurve(stars);
    renderer.render(scene, camera);
    
}
function animate2() {
    animationFrameId2=requestAnimationFrame(animate2);
    updateGalaxyGPU(stars, deltaTime, G, densityNeg, Radius_SImulation);
    controls.update();
    renderer.render(scene, camera);
    
}

function animateFree() {
    if (animationFrameId3) {
        cancelAnimationFrame(animationFrameId3);
    }
    if (animationFrameId2) {
        cancelAnimationFrame(animationFrameId2);
    }
    if (animationFrameId4) {
        cancelAnimationFrame(animationFrameId4);
    }
   
    animationFrameId= requestAnimationFrame(animateFree);
    controls.update();
    cloud.updateCloud();
    renderer.render(scene, camera);

    
}
let ReplaySimu = [];
let stepNow=0


let isFinished= false
function RunSimuWorker(){
    
    //console.log('totalstep',timestep)

    if (animationFrameId3) {
        cancelAnimationFrame(animationFrameId3);
    }
    if (animationFrameId) {
        cancelAnimationFrame(animationFrameId);
    }
    if (animationFrameId4) {
        cancelAnimationFrame(animationFrameId4);
    }


    animationFrameId2= requestAnimationFrame(RunSimuWorker);


    if (stepNow < timestep){
       // console.log("stepNow",stepNow)
        updateGalaxyGPULPReplay(stars, deltaTime, G, densityNeg, Radius_SImulation)
    } else{
        if(isFinished== false){
        //console.log("FINISHED ")
        isFinished = true
    }
        
    
    }

    stepNow++

    controls.update();
    renderer.render(scene, camera);


}

let timeNow=0

function Replay() {
    if (computeForcesKernel) {
      computeForcesKernel.destroy();
    }
    if (computeForcesKernel2) {
      computeForcesKernel2.destroy();
    }
    if (animationFrameId) {
      cancelAnimationFrame(animationFrameId);
    }
    if (animationFrameId2) {
      cancelAnimationFrame(animationFrameId2);
    }
    animationFrameId2 = requestAnimationFrame(Replay);
  
    if (timeNow < timestep) {
      const positions = new Float32Array(stars.length * 3);
      for (let j = 0; j < stars.length; j++) {
        stars[j].posX = ReplayArray[timeNow][j][0];
        stars[j].posY = ReplayArray[timeNow][j][1];
        stars[j].posZ = ReplayArray[timeNow][j][2];
        positions[j * 3] = stars[j].posX;
        positions[j * 3 + 1] = stars[j].posY;
        positions[j * 3 + 2] = stars[j].posZ;
      }
      cloud.geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
      timeNow++;
    } else {
      timeNow = 0;
    }
    controls.update();
    renderer.render(scene, camera);
  }
  

