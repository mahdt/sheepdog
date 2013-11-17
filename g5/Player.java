package sheepdog.g5;

import sheepdog.sim.Point;
import java.util.Arrays;

public class Player extends sheepdog.sim.Player {
    private int nblacks;
    private boolean mode;
    private int phase;
	private boolean dirFlags[];
	private int mult;
	private double speed;
	private Point prevDesired;
	private double constDist = 2;
	private double constDistX = 0.35, constDistY = 4;
	private final int up = 0, down = 1, left = 2, right = 3;
	

    public void init(int nblacks, boolean mode) {
        this.nblacks = nblacks;
        this.mode = mode;
        this.phase = 0;
        this.speed = 1;
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
    Point current = dogs[id-1];
	Point desired = new Point();
	//Point desired = new Point(99.5, 51);
	// phase 0 - dog in LH
	
	if (mode) {
		// advanced task code comes here
	} else {
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
	}

        return current;
    }

}