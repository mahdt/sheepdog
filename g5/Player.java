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
	private Point[] travelShape, undeliveredSheep;
	private LinkedList<Point> travelTraj;
	private boolean dogOnTraj, goCW, recompute;
	private Point desired, current, center;
	private int phase;
	private int tickCount, desiredTicks;


	// initialize our player
	public void init(int nblacks, boolean mode) {
		this.nblacks = nblacks;
		this.mode = mode;
		this.phase = 0;
		this.speed = 1.5;
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
		undeliveredSheep = computeUndeliveredSheep(sheeps);

		//Point desired = new Point(99.5, 51);
		// phase 0 - dog in LH
		
		// if the dog is moving to a new point
		if (dogOnTraj) {
			// if the desired point is reached
			if (current.equals(desired)) {
				dogOnTraj = false;
				if (phase == 0) phase++;
				if (phase == 1 && tickCount < desiredTicks) {
					dogOnTraj = true;
					tickCount++;
					
					// TODO: have dogs push inward at each vertex then return
					// to point and continue, instead of standing still
					return current;
				}
				if (phase == 1 && tickCount == desiredTicks) {
					tickCount = 0;
					dogOnTraj = false;
				} 
			}
			else {
				return (moveDog(current, desired));
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
		current = moveDog(current, desired);	

		/*if (mode) {
		// advanced task code comes here
	} else { }
		// basic task code comes here
		if (dogs.length == 1) {
			//if (phase != 4) phase = 3;
			if (phase == 0) {
				desired.x = 52;
				desired.y = 50;

				if (current.equals(desired)) {
					phase = 1;
					dirFlags[up] = true;
					System.out.println("entering phase 1");
				}
			}
			if (phase == 1 || phase == 2) {
				double offset = mult == 0 ? 0.5 : 0.0;
				if (dirFlags[up]) {				
					desired = new Point(current.x, constDist*mult + offset);
					if (current.equals(desired)) {
						dirFlags[up] = false;
						dirFlags[right] = true;
					}
				}
				else if (dirFlags[down]) {
					desired = new Point(current.x, 100.0 - constDist*mult - offset);
					if (current.equals(desired)) {
						dirFlags[down] = false;
						dirFlags[left] = true;
					}
				}
				else if (dirFlags[right]) {
					desired = new Point(100.0 - 0.5, current.y);
					if (current.equals(desired)) {
						dirFlags[right] = false;
						dirFlags[down] = true;
					}
				}
				else if (dirFlags[left]) {
					desired = new Point(50.0 + 0.5, current.y);
					if (current.equals(desired)) {
						dirFlags[left] = false;
						dirFlags[up] = true;
						if (phase ==1)
							mult++;
					}
				}
				if (phase == 2) {
					double minY = Double.MAX_VALUE;
					double maxY = Double.MIN_VALUE;
					for (int i = 0; i < sheeps.length; i++) {
						if (sheeps[i].y < minY && sheeps[i].y > 40.0 && sheeps[i].y < 60.0 && sheeps[i].x>50) minY = sheeps[i].y;
						if (sheeps[i].y > maxY && sheeps[i].y > 40.0 && sheeps[i].y < 60.0 && sheeps[i].x>50) maxY = sheeps[i].y;
					}
					double eps = 2.5;
					System.out.printf("maxy: %f, miny: %f\n", maxY, minY);
					if (maxY - minY <= eps && dirFlags[down]) {
						phase = 3;
						//speed = 0.5;
						//dirFlags[down] = false;
						System.out.println("Entering phase 3");
					}
				}
				if (mult * constDist >= 40.0)
					phase = 2;
			}
			if (phase == 3) {
				desired = new Point(99.5, 51.0);
				prevDesired = desired;
				Arrays.fill(dirFlags, false);
				if (current.equals(desired)) {
					phase = 4;
					mult = 0;
					dirFlags[up] = true;
				}
			}
			if (phase == 4) {
				if (dirFlags[up]) {
					desired = new Point(prevDesired.x - constDistX, prevDesired.y - constDistY);
					if (current.equals(desired)) {
						dirFlags[up] = false;
						dirFlags[down] = true; 
						prevDesired = desired;
					}
				}
				else if (dirFlags[down]) {
					desired = new Point(prevDesired.x - constDistX, prevDesired.y + constDistY);
					if (current.equals(desired)) {
						dirFlags[down] = false;
						dirFlags[up] = true; 
						prevDesired = desired;
					}
				}
				if (current.x < 51.0) phase = 5;
			}
			if (phase ==5)
				return current;
			current = moveDog (current, desired);
			System.out.println("desired " + desired.x + " " + desired.y + " up = "+dirFlags[up] + " and down = " + dirFlags[down]);
			System.out.println("current: " + current.x + " " + current.y);	
		}
	}*/

		return current;
	}
	
	// returns the sheep that have not been delivered
	public Point[] computeUndeliveredSheep(Point[] sheep) {
		ArrayList<Point> temp = new ArrayList<Point>();
		for (int i = 0; i < sheep.length; i++) {
			if (sheep[i].x > 50.0 + EPSILON) temp.add(sheep[i]);
		}
		return temp.toArray(new Point[temp.size()]);
	}
	
	
	// given the current point and desired point, computes the direction vector and
	// moves the dog towards that direction at the specified speed
	public Point moveDog(Point current, Point desired) {
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
