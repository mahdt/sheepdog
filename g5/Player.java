package sheepdog.g5;

import sheepdog.sim.Point;
import java.util.*;

public class Player extends sheepdog.sim.Player {

	//private boolean dirFlags[];
	//private int mult;
	//private Point prevDesired;
	//private double constDist = 2;
	//private double constDistX = 0.35, constDistY = 4;
	//private final int up = 0, down = 1, left = 2, right = 3;

	private final double EPSILON = 0.0001;
	private int nblacks;
	private boolean mode;
	private double speed;
	private Point[] prevDogs, travelShape, undeliveredSheep;
	private LinkedList<Point> travelTraj;
	private boolean dogOnTraj, goCW, recompute;
	private Point desired, current, center, lastVertex;
	private int phase;
	private int tickCount, desiredTicks;
	private int[] directions;
	
	public static double[][] allDogsInfo;


	// initialize our player
	public void init(int nblacks, boolean mode) {
		this.nblacks = nblacks;
		this.mode = mode;
		this.phase = 0;
		this.speed = 2;
		this.dogOnTraj = false;
		this.phase = 0;
		this.goCW = this.id % 2 == 0;
		this.recompute = true;
		this.tickCount = 0;
		this.desiredTicks = 3;
		this.center = new Point(50.0, 50.0);
		travelTraj = new LinkedList<Point>();
		//mult = 0;
		//dirFlags = new boolean[4];
		//Arrays.fill(dirFlags, false);
	}


	// returns the next position for the dog to move to
	public Point move(Point[] dogs, Point[] sheeps) {
		current = dogs[id - 1]; // our current position
		System.out.println("id:"+id);
		System.out.println("current.x:"+current.x);
		System.out.println("current.y:"+current.y);
		System.out.println("dog is in phase "+ phase);
		if (directions == null) directions = new int[dogs.length];
		
		undeliveredSheep = computeUndeliveredSheep(sheeps);

		//Point desired = new Point(99.5, 51);
		// phase 0 - dog in LH
		
		// if the dog is moving to a new point
		if (dogOnTraj) {
			// if the desired point is reached
			if (current.equals(desired)) {
				dogOnTraj = false;
				if (phase == 0) phase++;
				/*if (phase == 1 && tickCount < desiredTicks) {
					dogOnTraj = true;
					tickCount++;
					
					// TODO: have dogs push inward at each vertex then return
					// to point and continue, instead of standing still
					return current;
				}*/
				if (phase == 1) {
					lastVertex = current;
					phase++;
					speed = 0.5;
					desired = center;
					dogOnTraj = true;
				}
				else if (phase == 2) {
					calcDirections (dogs);
					//System.out.println("Should dog " + id + " wait: " + shouldDogWait());
					if (shouldDogWait(dogs, directions)) {
						dogOnTraj = true;
						return current;
					}
					phase--;
					speed = 2;
				}
				/*if (phase == 1 && tickCount == desiredTicks) {
					tickCount = 0;
					dogOnTraj = false;
				}*/
			}
			else if (phase == 2) {
				desired = lastVertex;
			}
			else {
				return (moveDog(current, desired, speed));
			}
		}

		// phase 0: move dog to gate
		if (phase == 0) {
			desired = center;
			dogOnTraj = true;
			recompute = true;
		}
		// phase 1, compute hull and begin traveling
		if (phase == 1) {
			// recomputes the hull for first traversal and each time a hull has been fully traversed
			if (recompute) {
				// compute the hull and grow it to get the dog's path
				travelShape = growHull(computeHull(undeliveredSheep));
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

		// move the dog to the next computed point
		current = moveDog(current, desired, speed);
		prevDogs = dogs.clone();

		return current;
	}
	
	public void calcDirections(Point[] dogs) {
		for (int i = 0; i < dogs.length; i++) {
			if (dogs[i].x < 50) directions[i] = -1;
			if ((dogs[i].y - prevDogs[i].y) > 0 || ((dogs[i].x - prevDogs[i].x) > 0 && dogs[i].y < 50) || 
													((dogs[i].x - prevDogs[i].x) < 0 && dogs[i].y > 50))
				directions[i] = 1;
			else if ((dogs[i].y - prevDogs[i].y) < 0 || ((dogs[i].x - prevDogs[i].x) > 0 && dogs[i].y > 50) || 
													((dogs[i].x - prevDogs[i].x) < 0 && dogs[i].y < 50))
				directions[i] = 0;			
		}
	}
	
	// returns the sheep that have not been delivered
	public Point[] computeUndeliveredSheep(Point[] sheep) {
		ArrayList<Point> temp = new ArrayList<Point>();
		for (int i = 0; i < sheep.length; i++) {
			if (sheep[i].x > 50.0 + EPSILON) temp.add(sheep[i]);
		}
		return temp.toArray(new Point[temp.size()]);
	}
	
	public boolean shouldDogWait(Point[] dogs, int[] directions) {
		System.out.println("In shouldDogWait()");
		int nextDogIndex = -1;
		if (goCW) {
			double closestYBelow = 0.0;
			for (int i = 0; i < dogs.length; i++) {
				if (directions[i] == 1 && dogs[i].y < current.y && dogs[i].y >= closestYBelow) {
					closestYBelow = dogs[i].y;
					nextDogIndex = i;
				}
			}
		}
		else {
			double closestYAbove = 0.0;
			for (int i = 0; i < dogs.length; i++) {
				System.out.println("i:" + i);
				if (directions[i] == 1 && dogs[i].y > current.y && dogs[i].y <= closestYAbove) {
					closestYAbove = dogs[i].y;
					nextDogIndex = i;
				}
			}
		}
		
		if (nextDogIndex == -1) return false;
		
		double d = computeHullDistance(current, new Point(dogs[nextDogIndex].x, dogs[nextDogIndex].y), goCW, travelShape);
		double hullLength = computeHullDistance(travelShape[0], travelShape[travelShape.length - 1], false, travelShape);
		System.out.println("d of dog " + id +":" + d);
		System.out.println("hullLength of dog " + id  + ":" + hullLength);
		System.out.println("maximum d for wait for dog " + id + ":" + (((double) dogs.length) / (2.0 * hullLength)));
			
		return (d < ((double) dogs.length) / (2.0 * hullLength));
	}
	
	// estimates the distance along a hull between two dogs
	// hull is defined as the hull of dog1.
	public double computeHullDistance(Point dog1, Point dog2, boolean direction, Point[] hull) {
		int startIndex = 0;
		for (int i = 0; i < hull.length; i++) {
			if (approxEqual(dog1, hull[i])) {
				startIndex = i;
				break;
			}
		}
		if (direction == false) { // counterclockwise
			double totalDist = 0;
			int i = 0;
			for (i = startIndex + 1; hull[i].y < dog2.y; i++) {
				totalDist += distance(hull[i - 1], hull[i]);
			}
			totalDist += distance(hull[i - 1], dog2);
			return totalDist;
		}
		else { // clockwise
			double totalDist = 0;
			int i = 0;
			for (i = startIndex - 1; hull[i].y > dog2.y; i--) {
				totalDist += distance(hull[i + 1], hull[i]);
			}
			totalDist += distance(hull[i + 1], dog2);
			return totalDist;
		}
	}
	
	public boolean approxEqual(Point a, Point b) {
		double eps = 0.0001;
		return Math.abs(a.x - b.x) <= eps && Math.abs(a.y - b.y) <= eps;
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
			return new Point(current.x + vector.x, current.y + vector.y);
		}
		// if dog is within 2 meters of desired point, move exactly to that point
		return desired;
	}
	
	public double distance(Point a, Point b) {
		return Math.sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
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
		minY = minY > 42.0 ? 42.0 : minY;
		maxY = maxY < 58.0 ? 58.0 : maxY;
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
			if(isCClockWise(stack.get(1).sheep, stack.get(0).sheep, nodes[i].sheep)) {
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
		Point[] newPoints = new Point[hull.length * 6];
		for (int i=0; i < hull.length; i++) {
			// original point
			newPoints[6 * i + 0] = new Point(hull[i]);
			// point above
			newPoints[6 * i + 1] = new Point(hull[i].x, hull[i].y - rad <= 0.0 ? 0.0 + EPSILON : hull[i].y - rad);
			// point below
			newPoints[6 * i + 2] = new Point(hull[i].x, hull[i].y + rad >= 100.0 ? 100.0 - EPSILON : hull[i].y + rad);
			// point on right
			newPoints[6 * i + 3] = new Point(hull[i].x + rad >= 100.0 ? 100.0 - EPSILON : hull[i].x + rad, hull[i].y);
			// top right
			newPoints[6 * i + 4] = new Point(hull[i].x + rad >= 100.0 ? 100.0 - EPSILON : hull[i].x + rad,
					hull[i].y - rad <= 0.0 ? 0.0 + EPSILON : hull[i].y - rad);
			// bottom right
			newPoints[6 * i + 5] = new Point(hull[i].x + rad >= 100.0 ? 100.0 - EPSILON : hull[i].x + rad,
					hull[i].y + rad >= 100.0 ? 100.0 - EPSILON : hull[i].y + rad);
		}
		return computeHull(newPoints);
	}

}
