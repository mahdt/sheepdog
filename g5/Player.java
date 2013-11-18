package sheepdog.g5;

import sheepdog.sim.Point;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.ArrayList;

public class Player extends sheepdog.sim.Player {
    private int nblacks;
    private boolean mode;
	private boolean dirFlags[];
	private int mult;
	private double speed;
	private Point prevDesired;
	private double constDist = 2;
	private double constDistX = 0.35, constDistY = 4;
	private final int up = 0, down = 1, left = 2, right = 3;
	private Point[] travelShape;
	private LinkedList<Point> travelTraj;
	private boolean dogOnTraj, goCW, recompute;
	private Point desired, current;
	private int phase;
	private int tickCount;
	
	class sheepNode implements Comparable<sheepNode>{
		public double distance;
		public double angle;
		public Point sheep;
		
		public sheepNode(Point P) {
			sheep = P;
		}
		
		public void setDist(double d) {
			distance =d;
		}
		
		public void setAngle(double a) {
			angle = a;
		}
		
		public void setDistAndAngle(Point p) {
			angle = Math.toDegrees(Math.atan2(p.y - sheep.y, sheep.x - p.x));
			distance = Math.sqrt((p.x - sheep.x)*(p.x - sheep.x) + (p.y - sheep.y)*(p.y - sheep.y)); 
		}
		
		public int compareTo (sheepNode other) {
			if (this.angle > other.angle) return 1;
			else if (this.angle < other.angle) return -1;
			else if (this.distance > other.distance) return 1;
			else if (this.distance < other.distance) return -1;
			else return 0;
		}
	}
	
	public boolean isCClockWise(Point p1, Point p2, Point p3) {
		double ux = (p2.x - p1.x);
		double uy = -(p2.y - p1.y);
		double vx = (p3.x - p1.x);
		double vy = -(p3.y - p1.y);
		return (ux * vy - uy * vx > 0);
	}
	
	public Point[] computeHull(Point[] sheep) {
		sheepNode[] nodes = new sheepNode[sheep.length + 2];
		double maxY = Double.MIN_VALUE;
		double minY = Double.MAX_VALUE;
		for (int i = 0; i < sheep.length; i++) {
			if (sheep[i].y < minY) minY = sheep[i].y;
			if (sheep[i].y > maxY) maxY = sheep[i].y;
		}
		minY = minY > 45.0 ? 45.0 : minY;
		maxY = maxY < 55.0 ? 55.0 : maxY;
		nodes[0] = new sheepNode(new Point(50.001, maxY));
		nodes[1] = new sheepNode(new Point(50.001, minY));
		
		sheepNode p0 = nodes[0];
		for (int i = 2; i<nodes.length; i++) {
			sheepNode temp = new sheepNode(sheep[i-2]);
			if (temp.sheep.y > p0.sheep.y || (temp.sheep.y == p0.sheep.y && temp.sheep.x > p0.sheep.x)) {
				p0 = temp;
			}
			nodes[i] = temp;
		}
		p0.setAngle(0.0);
		p0.setDist(0.0);
		
		for (int i = 0; i < nodes.length; i++) {
			nodes[i].setDistAndAngle(p0.sheep);
		}
		Arrays.sort(nodes);
		
		LinkedList<sheepNode> stack = new LinkedList<sheepNode>();
		stack.push(nodes[nodes.length-1]);
		stack.push(p0);
		int i = 1;
		while(i < nodes.length) {
			// check for clockwise turns since x-axis is flipped
			if(isCClockWise(stack.get(1).sheep, stack.get(0).sheep, nodes[i].sheep)) {
				stack.push(nodes[i]); 
				i++;
			}
			else {
				stack.pop();
			}
		}
		stack.pop();
		Point[] hull = new Point[stack.size()];
		for (int j=0; j<stack.size(); j++) {
			hull[j]=stack.get(j).sheep;
			System.out.println(hull[j].x + " "+ hull[j].y);
		}
		return hull;
	}
	
	public Point[] growHull(Point[] hull) {
		double rad = 2.0;
		Point[] newPoints = new Point[hull.length * 2];
		for (int i=0; i < hull.length; i++) {
			newPoints[2 * i] = new Point(hull[i].x + rad >= 100 ? 99.999 : hull[i].x + rad, hull[i].y - rad <= 0 ? 0.001 : hull[i].y - rad);
			newPoints[2 * i + 1] = new Point(hull[i].x + rad >= 100 ? 99.999 : hull[i].x + rad, hull[i].y + rad >= 100 ? 99.999 : hull[i].y + rad);
		}
		return computeHull(newPoints);
	}

    public void init(int nblacks, boolean mode) {
        this.nblacks = nblacks;
        this.mode = mode;
        this.phase = 0;
        this.speed = 1;
        this.dogOnTraj = false;
        this.phase = 0;
        this.goCW = false;
        this.recompute = true;
        this.tickCount = 0;
        travelTraj = new LinkedList<Point>();
        mult = 0;
        dirFlags = new boolean[4];
        Arrays.fill(dirFlags, false);
    }
	
	public Point moveDog(Point current, Point desired) {
		Point vector = new Point(desired.x - current.x, desired.y - current.y);
		double length = Math.sqrt(vector.x*vector.x + vector.y*vector.y);
		if (length == 0)
			return current;
		if (length > 2.0) {
			vector.x /= length/speed;
    		vector.y /= length/speed;
    		return new Point(current.x + vector.x, current.y + vector.y);
    	}
    	return desired;
	}
	

    // Return: the next position
    // my position: dogs[id-1]
    public Point move(Point[] dogs, // positions of dogs
                      Point[] sheeps) { // positions of the sheeps
    current = dogs[id-1];
	
	//Point desired = new Point(99.5, 51);
	// phase 0 - dog in LH
	if (dogOnTraj) {
		if (current.equals(desired)) {
			dogOnTraj = false;
			if (phase == 0) phase++;
			if (phase == 1 && tickCount < 3) {
				dogOnTraj = true;
				tickCount++;
				return current;
			}
			if (tickCount == 3) {
				tickCount = 0;
				dogOnTraj = false;
			} 
		}
		else
			return (moveDog(current, desired));
	}
	
	if (phase == 0) {
		desired = new Point(50,50);
		dogOnTraj = true;
		recompute = true;
	}
	if (phase == 1) {
		if (recompute) {
			ArrayList<Point> hullSheeps = new ArrayList<Point>();
			for (int i = 0; i < sheeps.length; i++) {
				if (sheeps[i].x > 50) hullSheeps.add(sheeps[i]);
			}
			travelShape = growHull(computeHull(hullSheeps.toArray(new Point[hullSheeps.size()])));
			travelTraj.clear();
			for (int i=0; i<travelShape.length; i++) {
				if (!goCW)
					travelTraj.push(travelShape[i]);
				else
					travelTraj.push(travelShape[travelShape.length - i -1]);
			}
			recompute = false;
		}
		
		if (travelTraj.size() == 0) {
			recompute = true;
			goCW = !goCW;
		}
		else {
			desired = travelTraj.pop();	
		}
		dogOnTraj = true;
	}
	
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

}
