package sheepdog.g5;

import sheepdog.sim.Point;
import java.util.*;

public class Player extends sheepdog.sim.Player {

	// constants
	private final int MAX_INITIAL_WAIT_TICKS = 60;
	private final double EPSILON = 0.000001;
	private final double FENCE_CENTER = 50.0;
	private final double FENCE_LEFT = 0.0;
	private final double FENCE_RIGHT = 100.0;
	private final double FENCE_TOP = 0.0;
	private final double FENCE_BOTTOM = 100.0;
	private final Point CENTER = new Point(FENCE_CENTER, FENCE_CENTER);
	private final Point LEFT_WALL = new Point(FENCE_LEFT, FENCE_CENTER);
	private final Point RIGHT_WALL = new Point(FENCE_RIGHT, FENCE_CENTER);
	private final double DOG_MAX_SPEED = 2.0;
	
	// tweakable parameters
	private final double SHEEP_FOLLOW_DISTANCE = 1.0; // distance behind a sheep a dog should try to stay in 1-by-1 strategy
	private final double HULL_GROWTH_RADIUS = 1.5; // distance to grow out the hull by
	private final double HULL_GROWTH_ENDGAME_MULT = 2.5; // multiplier on growth distance when endgame state is reached
	private final double HULL_TRAVERSAL_SPEED = 2.0; // how quickly the dog should traverse the hull
	private final double VERTEX_PUSH_DISTANCE = 0.85; // how far the dog should push towards the CENTER at each vertex of the hull
	private final double MIN_PUSH_X_VAL = 52.0; // the minimum x value that a dog pushing towards the CENTER can push to
	private final double ENDGAME_SHEEP_X_THRESHOLD = 51.0; // the maximum x value of all sheep that determines when to enter endgame state
	private final double ENDGAME_SPEED = 0.5; // how quickly the dog should travel in endgame state (normal mode)
	private final double MIN_VERTICAL_GATE_DIST = 8.0; // minimum vertical distance from the gate when computing the hull
	private final double MIN_VERTICAL_GATE_DIST_ENDGAME = 4.0; // same as above, but for endgame
	
	// instance variables
	private int nblacks, phase, ndogs, currentDelay, initialDelay, sheepGrabbed, chosenSheep;
	private double speed;
	private boolean mode, dogOnTraj, goCW, recompute, directionSet, targetingSheep, advDone, endGame, herdSet;
	private Point desired, current, lastVertex, desiredSheep;
	private int[] remainingSheep;
	private Point[] travelShape, undeliveredSheep, sheepLocations;
	private LinkedList<Point> travelTraj;


	// initialize our player
	public void init(int nblacks, boolean mode) {
		this.nblacks = nblacks;
		this.mode = mode;
		phase = 0;
		speed = 2;
		dogOnTraj = false;
		phase = 0;
		chosenSheep = -1;
		goCW = this.id % 2 == 0;
		recompute = true;
		currentDelay = 0;
		initialDelay = 0;
		travelTraj = new LinkedList<Point>();
		directionSet = false;
		targetingSheep = false;
		sheepGrabbed = -1;
		advDone = false;
		herdSet = false;
	}


	// returns the next position for the dog to move to
	public Point move(Point[] dogs, Point[] sheeps) {

		ndogs = dogs.length;
		current = dogs[id - 1];

		if (mode) {
			
			System.out.println("Phase for dog " + id + ": " + phase);
			
			if (phase == 0) {
				targetingSheep = false;
				desired = CENTER;
				if (approxEqual(current, CENTER)) {
					phase = 1;
				}
				return moveDog(current, desired, DOG_MAX_SPEED);
			}
			if (phase == 1) {
				if (!targetingSheep || sheeps[chosenSheep].x < FENCE_CENTER) {
					targetingSheep = true;
					remainingSheep = computeSheepOnSide(sheeps, true, false);
					if (remainingSheep == null) {
						phase = 2;
						return moveDog(current, CENTER, DOG_MAX_SPEED);
					}
					System.out.println("remaining sheep for dog " + id + ": " + remainingSheep.length);
					do {
						chosenSheep = remainingSheep[(int) (Math.random() * remainingSheep.length)];
					} while (sheepBeingPursued(dogs, current, sheeps[chosenSheep]) && remainingSheep.length > dogs.length);
					System.out.println("Chosen sheep for dog " + id + ": " + chosenSheep);
					if (current.x < FENCE_CENTER) return moveDog(current, CENTER, DOG_MAX_SPEED);
				}
				if (sheepBeingPursued(dogs, current, sheeps[chosenSheep])) {
					do {
						chosenSheep = remainingSheep[(int) (Math.random() * remainingSheep.length)];
					} while (sheepBeingPursued(dogs, current, sheeps[chosenSheep]) && remainingSheep.length > dogs.length);
				}
				desiredSheep = sheeps[chosenSheep];
				desired = getDirectionOfLength(SHEEP_FOLLOW_DISTANCE, CENTER, desiredSheep);
				desired.x += desiredSheep.x;
				desired.y += desiredSheep.y;
				return moveDog(current, desired, DOG_MAX_SPEED);
			}
			if (phase == 2) {
				targetingSheep = false;
				desired = CENTER;
				if (approxEqual(current, CENTER)) {
					phase = 3;
				}
				return moveDog(current, desired, DOG_MAX_SPEED);
			}
			if (phase == 3) {
				if (!targetingSheep || sheeps[chosenSheep].x > FENCE_CENTER) {
					targetingSheep = true;
					remainingSheep = computeSheepOnSide(sheeps, false, true);
					if (remainingSheep == null) {
						phase = 0;
						return moveDog(current, CENTER, DOG_MAX_SPEED);
					}
					do {
						chosenSheep = remainingSheep[(int) (Math.random() * remainingSheep.length)];
					} while (sheepBeingPursued(dogs, current, sheeps[chosenSheep]) && remainingSheep.length > dogs.length);
					if (current.x > FENCE_CENTER) return moveDog(current, CENTER, DOG_MAX_SPEED);
				}
				if (sheepBeingPursued(dogs, current, sheeps[chosenSheep])) {
					do {
						chosenSheep = remainingSheep[(int) (Math.random() * remainingSheep.length)];
					} while (sheepBeingPursued(dogs, current, sheeps[chosenSheep]) && remainingSheep.length > dogs.length);
				}
				desiredSheep = sheeps[chosenSheep];
				desired = getDirectionOfLength(SHEEP_FOLLOW_DISTANCE, CENTER, desiredSheep);
				desired.x += desiredSheep.x;
				desired.y += desiredSheep.y;
				return moveDog(current, desired, DOG_MAX_SPEED);
			}
			
			/*
			
			System.out.println("Dog ID " + id + " is in phase: " + phase);
			if (id > nblacks) return current;
			if (advDone && !isLastDog(dogs)) return moveDog(current, current.x > 50.0 ? RIGHT_WALL : LEFT_WALL, 2.0);
			int nBlacksLeft = computeUndeliveredSheepIndices(sheeps, true).size();
			if (nBlacksLeft == 0 && advDone && approxEqual(current, RIGHT_WALL)) {
				phase = 0;
				sheepGrabbed = -1;
				chosenSheep = -1;
				targetingSheep = false;
				return moveDog(current, CENTER, 2.0);
			}
			if (phase == 0) {
				System.out.println("dog " + id + " is running phase 0");
				desired = new Point(CENTER);
				System.out.println("IN PHASE 0");
				if (current.equals(CENTER)) {
					phase = advDone ? 2 : 1;
					System.out.println("GOING INTO PHASE 2");
				}
				else return moveDog(current, CENTER, 2.0);
			}
			else if (phase == 1) {
				//int nBlacksLeft = computeUndeliveredSheep(sheeps, true).length;
				//LinkedList<Integer> undeliveredIndices = computeUndeliveredSheepIndices(sheeps, true);
				//int nBlacksLeft = undeliveredIndices.size();
				//int nBlacksLeft = computeUndeliveredSheepIndices(sheeps, true).size();
				if (chosenSheep < sheeps.length && (!targetingSheep || sheeps[chosenSheep].x <= 50.0)) {
					targetingSheep = true;
					sheepGrabbed++;
					chosenSheep = (id - 1) + ndogs * sheepGrabbed;
					while(chosenSheep<sheeps.length && sheeps[chosenSheep].x <= 50){
						sheepGrabbed++;
						chosenSheep = (id - 1) + ndogs * sheepGrabbed;
					}
				}
				System.out.println("chosen Sheep for dog " + id + ":" + chosenSheep);
				if (chosenSheep >= nblacks) {
					advDone = true;
					if (isLastDog(dogs)) phase = 0;
					System.out.println("Moving to RIGHT_WALL");
					return moveDog(current, current.x > 50.0 ? RIGHT_WALL : LEFT_WALL, 2.0);
				}
				desiredSheep = sheeps[chosenSheep];
				desired = getDirectionOfLength(1.0, CENTER, desiredSheep);
				desired.x += desiredSheep.x;
				desired.y += desiredSheep.y;
				return moveDog(current, desired, 2.0);
			}

			else if (phase == 2 || phase == 3) {
			
				int chosenSheep = phase == 2 ? nblacks : 0;
				while (phase == 2 ? sheeps[chosenSheep].x >= 50.0 : sheeps[chosenSheep].x <= 50.0) {
					chosenSheep++;
					if (phase == 2 ? chosenSheep >= sheeps.length : chosenSheep >= nblacks) {
						phase = phase == 2 ? 3 : 2;
						return moveDog(current, CENTER, 2.0);
					}
				}
				desiredSheep = sheeps[chosenSheep];
				desired = getDirectionOfLength(1.0, CENTER, desiredSheep);
				desired.x += desiredSheep.x;
				desired.y += desiredSheep.y;
				return moveDog(current, desired, 2.0);
				
			}
			
			*/
		}

		else {
			
			//double thresholdRatio = -0.053 * (double) sheeps.length + 22.33;
			//double thresholdRatio = -0.030 * (double) sheeps.length + 19.24;
			double thresholdRatio = 330.26 * Math.pow(sheeps.length, -0.654);
			//double actualRatio = (double) computeUndeliveredSheepIndices(sheeps, false).size() / (double) dogs.length;
			double actualRatio = (double) sheeps.length / (double) dogs.length;
			
			boolean initialHerd = false;
			if (!herdSet) {
				initialHerd = actualRatio >= thresholdRatio;
				herdSet = true;
			}
			boolean herd = actualRatio >= thresholdRatio;
			
			System.out.println("threshold: " + thresholdRatio + ", actual: " + actualRatio);
			
			if (dogs.length == 1 || !herd) {
			
				if (phase == 0) {
					targetingSheep = false;
					desired = CENTER;
					if (approxEqual(current, CENTER)) {
						phase = 1;
					}
					return moveDog(current, desired, DOG_MAX_SPEED);
				}
				if (phase == 1) {
					if (!targetingSheep || sheeps[chosenSheep].x < FENCE_CENTER) {
						targetingSheep = true;
						remainingSheep = computeUndeliveredSheepNormal(sheeps);
						if (remainingSheep == null) {
							phase = 2;
							return moveDog(current, CENTER, DOG_MAX_SPEED);
						}
						System.out.println("remaining sheep for dog " + id + ": " + remainingSheep.length);
						do {
							chosenSheep = remainingSheep[(int) (Math.random() * remainingSheep.length)];
						} while (sheepBeingPursued(dogs, current, sheeps[chosenSheep]) && remainingSheep.length > dogs.length);
						System.out.println("Chosen sheep for dog " + id + ": " + chosenSheep);
						if (current.x < FENCE_CENTER) return moveDog(current, CENTER, DOG_MAX_SPEED);
					}
					if (sheepBeingPursued(dogs, current, sheeps[chosenSheep])) {
						do {
							chosenSheep = remainingSheep[(int) (Math.random() * remainingSheep.length)];
						} while (sheepBeingPursued(dogs, current, sheeps[chosenSheep]) && remainingSheep.length > dogs.length);
					}
					desiredSheep = sheeps[chosenSheep];
					desired = getDirectionOfLength(SHEEP_FOLLOW_DISTANCE, CENTER, desiredSheep);
					desired.x += desiredSheep.x;
					desired.y += desiredSheep.y;
					return moveDog(current, desired, DOG_MAX_SPEED);
				}
			
			/*
				if (phase == 0) {
					if (dogs.length != 1 && herd != initialHerd) {
						phase += 1;
						return current;
					}
					System.out.println("dog " + id + " is running phase 0");
					desired = new Point(CENTER);
					if (current.equals(CENTER)) {
						phase += 1;
					}
					else return moveDog(current, CENTER, DOG_MAX_SPEED);
				}
				else if (phase == 1) {
					//int nSheepLeft = computeUndeliveredSheepIndices(sheeps, false).size();
					if (chosenSheep < sheeps.length && (!targetingSheep || sheeps[chosenSheep].x <= FENCE_CENTER)) {
						targetingSheep = true;
						sheepGrabbed++;
						chosenSheep = (id - 1) + ndogs * sheepGrabbed;
						while(chosenSheep<sheeps.length && sheeps[chosenSheep].x <= FENCE_CENTER){
							sheepGrabbed++;
							chosenSheep = (id - 1) + ndogs * sheepGrabbed;
						}

					}
					System.out.println("chosen Sheep for dog " + id + ":" + chosenSheep);
					if (chosenSheep >= sheeps.length) {
						System.out.println("Moving to RIGHT_WALL");
						return moveDog(current, current.x > FENCE_CENTER ? RIGHT_WALL : LEFT_WALL, DOG_MAX_SPEED);
					}
					desiredSheep = sheeps[chosenSheep];
					desired = getDirectionOfLength(SHEEP_FOLLOW_DISTANCE, CENTER, desiredSheep);
					desired.x += desiredSheep.x;
					desired.y += desiredSheep.y;
					return moveDog(current, desired, DOG_MAX_SPEED);
				}
			*/
			}

			if (dogs.length % 2 != 0 && id > dogs.length / 2 + 1 && !directionSet) {
				goCW = !goCW;
				directionSet = true;
			}

			current = dogs[id - 1]; // our current position
			undeliveredSheep = computeUndeliveredSheep(sheeps, false);
			if (dogs.length < 3) initialDelay = 0;
			else {
				int posFromCenter = 0;
				if (dogs.length % 2 == 0) posFromCenter = Math.abs(id - dogs.length/2) - (id > dogs.length/2 ? 1 : 0);
				else posFromCenter = Math.abs(id - dogs.length/2 - 1);
				initialDelay = (MAX_INITIAL_WAIT_TICKS / (dogs.length / 2)) * posFromCenter;
				System.out.println("initialDelay for dog "+ id +" is: "+initialDelay+ " and pos from CENTER is " + posFromCenter);
			}
			double maxX = 0;                
			for (int i = 0; i < undeliveredSheep.length; i++) {
				if (undeliveredSheep[i].x > maxX)
					maxX = undeliveredSheep[i].x;
			}
			endGame = (maxX < ENDGAME_SHEEP_X_THRESHOLD);
			//Point desired = new Point(99.5, 51);
			// phase 0 - dog in LH

			// if the dog is moving to a new point
			if (dogOnTraj) {
				// if the desired point is reached
				if (current.equals(desired)) {
					dogOnTraj = false;
					if (phase == 0) phase++;
					speed = endGame ? ENDGAME_SPEED : DOG_MAX_SPEED;
				}
				else return (moveDog(current, desired, speed));
			}

			// phase 0: move dog to gate
			if (phase == 0) {
				if (currentDelay < initialDelay) {
					currentDelay++;
					return current;
				} 
				desired = CENTER;
				dogOnTraj = true;
				recompute = true;
			}
			// phase 1, compute hull and begin traveling
			if (phase == 1) {
				// recomputes the hull for first traversal and each time a hull has been fully traversed
				if (recompute) {
					// compute the hull and grow it to get the dog's path
					travelShape = computeFullPath(growHull(computeHull(undeliveredSheep)), VERTEX_PUSH_DISTANCE);
					//travelShape = growHull(computeHull(undeliveredSheep));
					travelTraj.clear();
					for (int i = 0; i < travelShape.length; i++) {
						if (!goCW) travelTraj.push(travelShape[i]);
						else travelTraj.push(travelShape[travelShape.length - 1 - i]);
					}
					recompute = false;
				}
				// if hull has been traversed, prepare to recompute and change direction
				if (travelTraj.size() == 0) {
					recompute = true;
					goCW = !goCW;
				}
				// otherwise, pop the next point to navigate to
				else {
					desired = travelTraj.pop();        
				}
				dogOnTraj = true;
			}
			if (phase == 2) {

			}
		}        
		// move the dog to the next computed point
		current = moveDog(current, desired, speed);        

		return current;
	}

	// returns the sheep that have not been delivered
	public Point[] computeUndeliveredSheep(Point[] sheep, boolean black) {
		ArrayList<Point> temp = new ArrayList<Point>();
		for (int i = 0; i < (black ? nblacks : sheep.length); i++) {
			if (sheep[i].x > FENCE_CENTER) temp.add(sheep[i]);
		}
		return temp.toArray(new Point[temp.size()]);
	}
	
	// returns the number of sheep of a specified color on a specified side
	public int[] computeSheepOnSide(Point[] sheep, boolean black, boolean left) {
		ArrayList<Integer> temp = new ArrayList<Integer>();
		for (int i = black ? 0 : nblacks; i < (black ? nblacks : sheep.length); i++) {
			if (left ? sheep[i].x <= FENCE_CENTER : sheep[i].x >= FENCE_CENTER) temp.add(i);
		}
		if (temp.isEmpty()) return null;
		int[] returnVal = new int[temp.size()];
		for (int i = 0; i < returnVal.length; i++) returnVal[i] = temp.get(i).intValue();
		return returnVal;
	}

	public Point[] computeDeliveredSheepOfOneColor(Point[] sheep, boolean black) {
		ArrayList<Point> temp = new ArrayList<Point>();
		for (int i = black ? 0 : nblacks; i < (black ? nblacks : sheep.length); i++) {
			if (sheep[i].x <= FENCE_CENTER) temp.add(sheep[i]);
		}
		return temp.toArray(new Point[temp.size()]);
	}

	// returns the sheep that have not been delivered
	public LinkedList<Integer> computeUndeliveredSheepIndices(Point[] sheep, boolean black) {
		LinkedList<Integer> temp = new LinkedList<Integer>();
		for (int i = 0; i < (black ? nblacks : sheep.length); i++) {
			if (sheep[i].x > FENCE_CENTER) temp.add(i);
		}
		return temp;
	}
	
	// returns the indices of sheep that have not been delivered (used in normal mode)
	public int[] computeUndeliveredSheepNormal(Point[] sheep) {
		ArrayList<Integer> temp = new ArrayList<Integer>();
		for (int i = 0; i < sheep.length; i++) {
			if (sheep[i].x > FENCE_CENTER) temp.add(i);
		}
		if (temp.isEmpty()) return null;
		int[] returnVal = new int[temp.size()];
		for (int i = 0; i < returnVal.length; i++) returnVal[i] = temp.get(i).intValue();
		return returnVal;
	}

	public boolean sheepIsAssigned(int sheepIndex, Point[] sheeps, Point[] dogs) {
		for(int dogIndex = 0;dogIndex<dogs.length;dogIndex++) {
			if(distance(sheeps[sheepIndex],dogs[dogIndex]) < 2.0 && dogIndex+1!=id){
				return true;
			}
		}
		return false;
	}
	
	public boolean sheepBeingPursued(Point[] dogs, Point current, Point sheep) {
		for (int i = 0; i < dogs.length; i++) {
			if (approxEqual(getDirectionOfLength(1.0, CENTER, dogs[i]), getDirectionOfLength(1.0, CENTER, sheep))
					&& distance(dogs[i], sheep) < distance(current, sheep)) {
				return true;
			}
		}
		return false;
	}

	public boolean isLastDog(Point[] dogs){
		System.out.println("In isLastDog()");
		for(int i = 0;i<dogs.length;i++){
			System.out.println("In for loop isLastDog()");
			if(id!=i+1 && !approxEqual(dogs[i], RIGHT_WALL)){
				System.out.println("dog "+ i + " is not at the wall");
				return false;
			}
		}
		return true;
	}


	// given the current point and desired point, computes the direction vector and
	// moves the dog towards that direction at the specified speed
	public Point moveDog(Point current, Point desired, double speed) {
		Point vector = new Point(desired.x - current.x, desired.y - current.y);
		double length = Math.sqrt(vector.x*vector.x + vector.y*vector.y);
		if (length == 0) return current;
		if (length > 2.0) {
			vector.x /= length; vector.x *= speed - EPSILON;
			vector.y /= length; vector.y *= speed - EPSILON;
			return clamp(new Point(current.x + vector.x, current.y + vector.y), FENCE_TOP, FENCE_BOTTOM, FENCE_LEFT, FENCE_RIGHT);
		}
		// if dog is within 2 meters of desired point, move exactly to that point
		return clamp(desired, FENCE_TOP, FENCE_BOTTOM, FENCE_LEFT, FENCE_RIGHT);
	}

	public Point clamp(Point p, double up, double down, double left, double right) {
		Point n = new Point(p);
		if (n.x > right - EPSILON) n.x = right - EPSILON;
		else if (n.x < left + EPSILON) n.x = left + EPSILON;
		if (n.y < up + EPSILON) n.y = up + EPSILON;
		else if (n.y > down - EPSILON) n.y = down - EPSILON;
		return n;
	}

	// computes unit vector between two points
	public Point getDirectionOfLength(double length, Point current, Point desired) {
		Point vector = new Point(desired.x - current.x, desired.y - current.y);
		double mag = Math.sqrt(vector.x*vector.x + vector.y*vector.y);
		vector.x /= mag; vector.x *= length;
		vector.y /= mag; vector.y *= length;
		return new Point(vector.x, vector.y);
	}

	// distance
	public double distance(Point current, Point desired) {
		Point vector = new Point(desired.x - current.x, desired.y - current.y);
		return Math.sqrt(vector.x*vector.x + vector.y*vector.y);
	}

	public boolean approxEqual(Point a, Point b) {
		if (Math.abs(a.x - b.x) > EPSILON) return false;
		if (Math.abs(a.y - b.y) > EPSILON) return false;
		return true;
	}

	// inner class to represent a node in our convex hull computation
	class SheepNode implements Comparable<SheepNode>{
		public double distance, angle;
		public Point sheep;
		public SheepNode(Point p) { sheep = p; }
		public void setDist(double d) { distance = d; }
		public void setAngle(double a) { angle = a; }
		public void setDistAndAngle(Point p) {
			angle = Math.toDegrees(Math.atan2(p.y - sheep.y, sheep.x - p.x));
			distance = Math.sqrt((p.x - sheep.x)*(p.x - sheep.x) + (p.y - sheep.y)*(p.y - sheep.y));
		}
		public int compareTo (SheepNode other) {
			if (this.angle > other.angle) return 1;
			else if (this.angle < other.angle) return -1;
			else if (this.distance > other.distance) return 1;
			else if (this.distance < other.distance) return -1;
			else return 0;
		}
	}

	// determines if three given points form a counterclockwise turn
	public boolean isCClockWise(Point p1, Point p2, Point p3) {
		double ux = (p2.x - p1.x);
		double uy = -(p2.y - p1.y);
		double vx = (p3.x - p1.x);
		double vy = -(p3.y - p1.y);
		return (ux * vy - uy * vx > 0);
	}

	// uses Graham's algorithm to compute the convex hull given the coordinates of all undelivered sheep
	public Point[] computeHull(Point[] sheep) {
		SheepNode[] nodes = new SheepNode[sheep.length + 2];

		// compute max and min y values to determine fencepoint nodes
		double maxY = Double.MIN_VALUE;
		double minY = Double.MAX_VALUE;
		for (int i = 0; i < sheep.length; i++) {
			if (sheep[i].y < minY) minY = sheep[i].y;
			if (sheep[i].y > maxY) maxY = sheep[i].y;
		}
		// if all the sheep are above or below the gate, adjust the max/min values accordingly
		double distFromGate = endGame ? MIN_VERTICAL_GATE_DIST_ENDGAME : MIN_VERTICAL_GATE_DIST;
		minY = minY > FENCE_CENTER - distFromGate ? FENCE_CENTER - distFromGate : minY;
		maxY = maxY < FENCE_CENTER + distFromGate ? FENCE_CENTER + distFromGate : maxY;
		nodes[0] = new SheepNode(new Point(FENCE_CENTER, maxY)); // exceptions might happen here, put + EPSILON back if so
		nodes[1] = new SheepNode(new Point(FENCE_CENTER, minY));

		// find p0, the rightmost lowest point, and set distance and angle to zero
		SheepNode p0 = nodes[0];
		for (int i = 2; i<nodes.length; i++) {
			SheepNode temp = new SheepNode(sheep[i-2]);
			if (temp.sheep.y > p0.sheep.y || (temp.sheep.y == p0.sheep.y && temp.sheep.x > p0.sheep.x)) {
				p0 = temp;
			}
			nodes[i] = temp;
		}
		p0.setAngle(0.0);
		p0.setDist(0.0);

		// set distance and angle of all other points relative to p0, and sort them
		for (int i = 0; i < nodes.length; i++) {
			nodes[i].setDistAndAngle(p0.sheep);
		}
		Arrays.sort(nodes);

		// use Graham's to compute the convex hull. at the end of the loop, the stack will contain the hull points
		LinkedList<SheepNode> stack = new LinkedList<SheepNode>();
		stack.push(nodes[nodes.length-1]);
		stack.push(p0);
		int i = 1;
		while(i < nodes.length) {
			if (isCClockWise(stack.get(1).sheep, stack.get(0).sheep, nodes[i].sheep)) {
				stack.push(nodes[i]); 
				i++;
			}
			else {
				stack.pop();
			}
		}
		stack.pop(); // to avoid p0 on stack twice

		// convert to Point array
		Point[] hull = new Point[stack.size()];
		for (int j=0; j<stack.size(); j++) {
			hull[j]=stack.get(j).sheep;
			System.out.println(hull[j].x + " "+ hull[j].y); // for debugging
		}
		return hull;
	}

	// grows the hull outward so the dog will enclose all the sheep
	public Point[] growHull(Point[] hull) {
		double rad = HULL_GROWTH_RADIUS;
		double rightRad = endGame ? rad * HULL_GROWTH_ENDGAME_MULT : rad;
		Point[] newPoints = new Point[hull.length * 6];
		for (int i=0; i < hull.length; i++) {
			// original point
			newPoints[6 * i + 0] = new Point(hull[i]);
			// point above
			newPoints[6 * i + 1] = new Point(hull[i].x, hull[i].y - rad <= FENCE_TOP ? FENCE_TOP + EPSILON : hull[i].y - rad);
			// point below
			newPoints[6 * i + 2] = new Point(hull[i].x, hull[i].y + rad >= FENCE_BOTTOM ? FENCE_BOTTOM - EPSILON : hull[i].y + rad);
			// point on right
			newPoints[6 * i + 3] = new Point(hull[i].x + rightRad >= FENCE_RIGHT ? FENCE_RIGHT - EPSILON : hull[i].x + rightRad, hull[i].y);
			// top right
			newPoints[6 * i + 4] = new Point(hull[i].x + rightRad >= FENCE_RIGHT ? FENCE_RIGHT - EPSILON : hull[i].x + rightRad,
					hull[i].y - rad <= FENCE_TOP ? FENCE_TOP + EPSILON : hull[i].y - rad);
			// bottom right
			newPoints[6 * i + 5] = new Point(hull[i].x + rightRad >= FENCE_RIGHT ? FENCE_RIGHT - EPSILON : hull[i].x + rightRad,
					hull[i].y + rad >= FENCE_BOTTOM ? FENCE_BOTTOM - EPSILON : hull[i].y + rad);
		}
		return computeHull(newPoints);
	}

	public Point[] computeFullPath(Point[] grownHull, double moveInDistance) {
		ArrayList<Point> fullPath = new ArrayList<Point>();
		for (int i = 0; i < grownHull.length; i++) {
			fullPath.add(grownHull[i]);
			Point inPoint = getDirectionOfLength(moveInDistance, grownHull[i], CENTER);
			inPoint.x += grownHull[i].x;
			inPoint.y += grownHull[i].y;
			fullPath.add(clamp(inPoint, FENCE_TOP, FENCE_BOTTOM, MIN_PUSH_X_VAL, FENCE_RIGHT));
			fullPath.add(grownHull[i]);
		}
		return fullPath.toArray(new Point[fullPath.size()]);
	}

}
