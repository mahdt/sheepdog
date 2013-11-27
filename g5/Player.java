package sheepdog.g5;

import sheepdog.sim.Point;
import java.util.*;

public class Player extends sheepdog.sim.Player {

	private final double EPSILON = 0.0001;
	private int nblacks;
	private boolean mode;
	private double speed;
	private Point[] travelShape, undeliveredSheep;
	private LinkedList<Point> travelTraj;
	private boolean dogOnTraj, goCW, recompute, directionSet;
	private Point desired, current, center, lastVertex, desiredSheep;
	private int phase;
	private int currentDelay, initialDelay, chosenSheep;
	private static final int maxTicks = 60;
	private boolean endGame;


	// initialize our player
	public void init(int nblacks, boolean mode) {
		this.nblacks = nblacks;
		this.mode = mode;
		this.phase = 0;
		this.speed = 2;
		this.dogOnTraj = false;
		this.phase = 0;
		this.chosenSheep = -1;
		this.goCW = this.id % 2 == 0;
		this.recompute = true;
		this.currentDelay = 0;
		this.initialDelay = 0;
		this.center = new Point(50.0, 50.0);
		travelTraj = new LinkedList<Point>();
		this.directionSet = false;
		//mult = 0;
		//dirFlags = new boolean[4];
		//Arrays.fill(dirFlags, false);
	}


	// returns the next position for the dog to move to
	public Point move(Point[] dogs, Point[] sheeps) {
		
		current = dogs[id - 1];
		// advanced mode?
		if (mode) {
		
			if (phase == 0) {
				desired = new Point(center);
				if (current.equals(center)) phase++;
				else return moveDog(current, center, 2.0 - EPSILON);
			}
			
			else {
				if (chosenSheep > -1) {
					if (getSide(sheeps[chosenSheep].x, sheeps[chosenSheep].y) == 0) chosenSheep = -1;
				}
			
				if (chosenSheep < 0) {
					undeliveredSheep = computeUndeliveredSheep(sheeps, true);
					chosenSheep = (int) (Math.random() * undeliveredSheep.length);				
				}
			
				desiredSheep = moveSheep(undeliveredSheep, chosenSheep);
			
				desired = getUnitDirection(center, desiredSheep);
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
					if (phase == 0) {
						phase++;
					}
					if (phase == 1) {
						lastVertex = current;
						phase++;
						speed = endGame ? 0.5 : 1.0;
						desired = center;
						dogOnTraj = true;
					}
					else if (phase == 2) {
						phase--;
						speed = endGame ? 0.75 : 2;
					}
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
		}	
		// move the dog to the next computed point
		current = moveDog(current, desired, speed);	

		return current;
	}
	
	// returns the sheep that have not been delivered
	public Point[] computeUndeliveredSheep(Point[] sheep, boolean black) {
		ArrayList<Point> temp = new ArrayList<Point>();
		for (int i = 0; i < (black ? nblacks : sheep.length); i++) {
			if (sheep[i].x > 50.0 + EPSILON) temp.add(sheep[i]);
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
			return new Point(current.x + vector.x, current.y + vector.y);
		}
		// if dog is within 2 meters of desired point, move exactly to that point
		return desired;
	}
	
	// computes unit vector between two points
	public Point getUnitDirection(Point current, Point desired) {
		Point vector = new Point(desired.x - current.x, desired.y - current.y);
		double length = Math.sqrt(vector.x*vector.x + vector.y*vector.y);
		vector.x /= length;
		vector.y /= length;
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
	
	
	
	
	Point moveSheep(Point[] sheeps, int sheepId) {
        Point thisSheep = sheeps[sheepId];
        double dspeed = 0;
        System.out.println("current: ");
        System.out.println(current.x + ", " + current.y);
        Point closestDog = current;
        double dist = distance(thisSheep, closestDog);
        assert dist > 0;

        if (dist < 2.0) // run dist
            dspeed = 1.0; // run speed
        else if (dist < 10.0) // walk dist
            dspeed = 0.1; // walk speed
        
        // offset from dogs
        double ox_dog = (thisSheep.x - closestDog.x) / dist * dspeed;
        double oy_dog = (thisSheep.y - closestDog.y) / dist * dspeed;

        // offset from clustering
        double ox_cluster = 0, oy_cluster = 0;

        // aggregate offsets then normalize
        for (int i = 0; i < sheeps.length; ++i) {
            // skip this sheep itself
            if (i == sheepId) continue;

            double d = distance(thisSheep, sheeps[i]);

            // ignore overlapping sheep
            if (d < 1.0 && d > 0) { // 1.0 = CLUSTER_DIST
                // add an unit vector to x-axis, y-axis
                ox_cluster += ((thisSheep.x - sheeps[i].x) / d);
                oy_cluster += ((thisSheep.y - sheeps[i].y) / d);
            }
        }
        // normalize by length
        double l = Math.sqrt(ox_cluster * ox_cluster + oy_cluster * oy_cluster);
        if (l != 0) {
            ox_cluster = ox_cluster / l * 0.05; // 0.05 = CLUSTER_SPEED
            oy_cluster = oy_cluster / l * 0.05;
        }

        double ox = ox_dog + ox_cluster, oy = oy_dog + oy_cluster;
        
        Point npos = updatePosition(thisSheep, ox, oy);

        return npos;

    }

    // update the current point according to the offsets
    Point updatePosition(Point now, double ox, double oy) {
        double nx = now.x + ox, ny = now.y + oy;
		double dimension = 100.0;
        // hit the left fence        
        if (nx < 0) {
            //            System.err.println("SHEEP HITS THE LEFT FENCE!!!");

            // move the point to the left fence
            Point temp = new Point(0, now.y);
            // how much we have already moved in x-axis?
            double moved = 0 - now.x;
            // how much we still need to move
            // BUT in opposite direction
            double ox2 = -(ox - moved); 
            return updatePosition(temp, ox2, oy);
        }
        // hit the right fence
        if (nx > dimension) {
            //            System.err.println("SHEEP HITS THE RIGHT FENCE!!!");

            // move the point to the right fence
            Point temp = new Point(dimension, now.y);
            double moved = (dimension - now.x);
            double ox2 = -(ox - moved);
            return updatePosition(temp, ox2, oy);
        }
        // hit the up fence
        if (ny < 0) {
            //            System.err.println("SHEEP HITS THE UP FENCE!!!");

            // move the point to the up fence
            Point temp = new Point(now.x, 0);
            double moved = 0 - now.y;
            double oy2 = -(oy - moved);
            return updatePosition(temp, ox, oy2);
        }
        // hit the bottom fence
        if (ny > dimension) {
            //            System.err.println("SHEEP HITS THE BOTTOM FENCE!!!");

            Point temp = new Point(now.x, dimension);
            double moved = (dimension - now.y);
            double oy2 = -(oy - moved);
            return updatePosition(temp, ox, oy2);
        }

        assert nx >= 0 && nx <= dimension;
        assert ny >= 0 && ny <= dimension;
        
        // hit the middle fence
        if (hitTheFence(now.x, now.y, nx, ny)) {
            //            System.err.println("SHEEP HITS THE CENTER FENCE!!!");
            //            System.err.println(nx + " " + ny);
            //            System.err.println(ox + " " + oy);

            // move the point to the fence
            Point temp = new Point(50, now.y);
            double moved = (50 - now.x);
            double ox2 = -(ox-moved);
            return updatePosition(temp, ox2, oy);
        }

        // otherwise, we are good
        return new Point(nx, ny);
    }
	
    boolean hitTheFence(double x1, double y1,
                        double x2, double y2) {
        double dimension = 100.0;
        // on the same side
        if (getSide(x1,y1) == getSide(x2, y2))
            return false;

        // one point is on the fence
        if (getSide(x1,y1) == 2 || getSide(x2,y2) == 2)
            return false;
        
        // compute the intersection with (50, y3)
        // (y3-y1)/(50-x1) = (y2-y1)/(x2-x1)

        double y3 = (y2-y1)/(x2-x1)*(50-x1)+y1;

        assert y3 >= 0 && y3 <= dimension;

        // pass the openning?
        if (y3 >= 49.0 && y3 <= 51.0)
            return false;
        else
            return true;
    }
	
    int getSide(double x, double y) {
   		double dimension = 100.0;
        if (x < dimension * 0.5)
            return 0;
        else if (x > dimension * 0.5)
            return 1;
        else
            return 2;
    }
	
	
	

}
