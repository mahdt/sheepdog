package sheepdog.g5;

import sheepdog.sim.Point;
import java.util.*;

public class Player extends sheepdog.sim.Player {

	private final double EPSILON = 0.000001;
	private int nblacks;
	private boolean mode;
	private double speed;
	private Point[] travelShape, undeliveredSheep;
	private LinkedList<Point> travelTraj;
	private boolean dogOnTraj, goCW, recompute, directionSet, targetingSheep, advDone;
	private Point desired, current, center, frontWall, backWall, lastVertex, desiredSheep;
	private int phase;
	private int ndogs, currentDelay, initialDelay, sheepGrabbed, chosenSheep;
	private static final int maxTicks = 60;
	private boolean endGame;


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
		center = new Point(50.0, 50.0);
		frontWall = new Point(0.0 + EPSILON, 50.0);
		backWall = new Point(100.0 - EPSILON, 50.0);
		travelTraj = new LinkedList<Point>();
		directionSet = false;
		targetingSheep = false;
		sheepGrabbed = -1;
		advDone = false;
	}


	// returns the next position for the dog to move to
	public Point move(Point[] dogs, Point[] sheeps) {
		
		ndogs = dogs.length;
		current = dogs[id - 1];
		
		// TODO: If there are more dogs than black sheep, move out the center dogs instead of assigning by ID
		// TODO: Perhaps task dogs who are "done" to try to push white sheep out of the way
		// TODO: If any white sheep make it across by accident, find a way to get them back
		// TODO: If a black sheep is inadvertently delivered by a dog who was not targeting it, exceptions will happen
		// TODO: Use both modes to empirically determine the dog/sheep ratio at which 1-by-1 is faster than convex hull
		if (mode) {
			
			if (id > nblacks) return current;
			if (advDone) return moveDog(current, current.x > 50.0 ? backWall : frontWall, 2.0 - EPSILON);
			
			if (phase == 0) {
				desired = new Point(center);
				if (current.equals(center)) phase++;
				else return moveDog(current, center, 2.0 - EPSILON);
			}
			
			else {
				int nBlacksLeft = computeUndeliveredSheep(sheeps, true).length;
				if (!targetingSheep || sheeps[chosenSheep].x <= 50.0) {
					targetingSheep = true;
					sheepGrabbed++;
					chosenSheep = (id - 1) + ndogs * sheepGrabbed;
					if (chosenSheep >= nblacks) {
						advDone = true;
						return moveDog(current, current.x > 50.0 ? backWall : frontWall, 2.0 - EPSILON);
					}
				}
	
				//desiredSheep = moveSheep(undeliveredSheep, chosenSheep);
				desiredSheep = sheeps[chosenSheep];
				desired = getDirectionOfLength(1.0, center, desiredSheep);
				desired.x += desiredSheep.x;
				desired.y += desiredSheep.y;
				return moveDog(current, desired, 2.0 - EPSILON);
			}
		}
		
		else {
		
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
				initialDelay = (maxTicks/(dogs.length/2))*(posFromCenter);
				System.out.println("initialDelay for dog "+ id +" is: "+initialDelay+ " and pos from center is " + posFromCenter);
			}
			double maxX = 0;		
			for (int i = 0; i < undeliveredSheep.length; i++) {
				if (undeliveredSheep[i].x > maxX)
					maxX = undeliveredSheep[i].x;
			}
			endGame = (maxX < 52);
			//Point desired = new Point(99.5, 51);
			// phase 0 - dog in LH
		
			// if the dog is moving to a new point
			if (dogOnTraj) {
				// if the desired point is reached
				if (current.equals(desired)) {
					dogOnTraj = false;
					if (phase == 0) phase++;
					speed = endGame ? 0.5 : 2;
				}
				else return (moveDog(current, desired, speed));
			}

			// phase 0: move dog to gate
			if (phase == 0) {
				if (currentDelay < initialDelay) {
					currentDelay++;
					return current;
				} 
				desired = center;
				dogOnTraj = true;
				recompute = true;
			}
			// phase 1, compute hull and begin traveling
			if (phase == 1) {
				// recomputes the hull for first traversal and each time a hull has been fully traversed
				if (recompute) {
					// compute the hull and grow it to get the dog's path
					travelShape = computeFullPath(growHull(computeHull(undeliveredSheep)), 0.75);
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
			if (sheep[i].x > 50.0) temp.add(sheep[i]);
		}
		return temp.toArray(new Point[temp.size()]);
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
			return clamp(new Point(current.x + vector.x, current.y + vector.y), 0.0, 100.0, 0.0, 100.0);
		}
		// if dog is within 2 meters of desired point, move exactly to that point
		return clamp(desired, 0.0, 100.0, 0.0, 100.0);
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
		double distFromGate = endGame ? 4.0 : 8;
		minY = minY > 50 - distFromGate ? 50 - distFromGate : minY;
		maxY = maxY < 50 + distFromGate ? 50 + distFromGate : maxY;
		nodes[0] = new SheepNode(new Point(50.0 + EPSILON, maxY));
		nodes[1] = new SheepNode(new Point(50.0 + EPSILON, minY));

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
		double rad = 1.5;
		double rightRad = endGame ? rad * 2.5 : rad;
		Point[] newPoints = new Point[hull.length * 6];
		for (int i=0; i < hull.length; i++) {
			// original point
			newPoints[6 * i + 0] = new Point(hull[i]);
			// point above
			newPoints[6 * i + 1] = new Point(hull[i].x, hull[i].y - rad <= 0.0 ? 0.0 + EPSILON : hull[i].y - rad);
			// point below
			newPoints[6 * i + 2] = new Point(hull[i].x, hull[i].y + rad >= 100.0 ? 100.0 - EPSILON : hull[i].y + rad);
			// point on right
			newPoints[6 * i + 3] = new Point(hull[i].x + rightRad >= 100.0 ? 100.0 - EPSILON : hull[i].x + rightRad, hull[i].y);
			// top right
			newPoints[6 * i + 4] = new Point(hull[i].x + rightRad >= 100.0 ? 100.0 - EPSILON : hull[i].x + rightRad,
					hull[i].y - rad <= 0.0 ? 0.0 + EPSILON : hull[i].y - rad);
			// bottom right
			newPoints[6 * i + 5] = new Point(hull[i].x + rightRad >= 100.0 ? 100.0 - EPSILON : hull[i].x + rightRad,
					hull[i].y + rad >= 100.0 ? 100.0 - EPSILON : hull[i].y + rad);
		}
		return computeHull(newPoints);
	}
	
	public Point[] computeFullPath(Point[] grownHull, double moveInDistance) {
		ArrayList<Point> fullPath = new ArrayList<Point>();
		
		for (int i = 0; i < grownHull.length; i++) {
			fullPath.add(grownHull[i]);
			Point inPoint = getDirectionOfLength(moveInDistance, grownHull[i], center);
			inPoint.x += grownHull[i].x;
			inPoint.y += grownHull[i].y;
			fullPath.add(clamp(inPoint, 0.0, 100.0, 52.0, 100.0));
			fullPath.add(grownHull[i]);
		}
		
		return fullPath.toArray(new Point[fullPath.size()]);
	}

}
