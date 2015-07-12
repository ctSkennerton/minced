
import java.util.ArrayList;
import java.util.List;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Map.Entry;

/**
 * The Node class contains the interval tree information for one single node
 *
 * @author Kevin Dolan
 */
public class IntervalNode<Type> {

	private SortedMap<Interval<Type>, List<Interval<Type>>> intervals;
	private long center;
	private IntervalNode<Type> leftNode;
	private IntervalNode<Type> rightNode;

	public IntervalNode() {
		intervals = new TreeMap<Interval<Type>, List<Interval<Type>>>();
		center = 0;
		leftNode = null;
		rightNode = null;
	}

	public IntervalNode(List<Interval<Type>> intervalList) {

		intervals = new TreeMap<Interval<Type>, List<Interval<Type>>>();

		SortedSet<Long> endpoints = new TreeSet<Long>();

		for(Interval<Type> interval: intervalList) {
			endpoints.add(interval.getStart());
			endpoints.add(interval.getEnd());
		}

		long median = getMedian(endpoints);
		center = median;

		List<Interval<Type>> left = new ArrayList<Interval<Type>>();
		List<Interval<Type>> right = new ArrayList<Interval<Type>>();

		for(Interval<Type> interval : intervalList) {
			if(interval.getEnd() < median)
				left.add(interval);
			else if(interval.getStart() > median)
				right.add(interval);
			else {
				List<Interval<Type>> posting = intervals.get(interval);
				if(posting == null) {
					posting = new ArrayList<Interval<Type>>();
					intervals.put(interval, posting);
				}
				posting.add(interval);
			}
		}

		if(left.size() > 0)
			leftNode = new IntervalNode<Type>(left);
		if(right.size() > 0)
			rightNode = new IntervalNode<Type>(right);
	}

	/**
	 * Perform a stabbing query on the node
	 * @param time the time to query at
	 * @return	   all intervals containing time
	 */
	public List<Interval<Type>> stab(long time) {
		List<Interval<Type>> result = new ArrayList<Interval<Type>>();

		for(Entry<Interval<Type>, List<Interval<Type>>> entry : intervals.entrySet()) {
			if(entry.getKey().contains(time))
				for(Interval<Type> interval : entry.getValue())
					result.add(interval);
			else if(entry.getKey().getStart() > time)
				break;
		}

		if(time < center && leftNode != null)
			result.addAll(leftNode.stab(time));
		else if(time > center && rightNode != null)
			result.addAll(rightNode.stab(time));
		return result;
	}

	/**
	 * Perform an interval intersection query on the node
	 * @param target the interval to intersect
	 * @return		   all intervals containing time
	 */
	public List<Interval<Type>> query(Interval<?> target) {
		List<Interval<Type>> result = new ArrayList<Interval<Type>>();

		for(Entry<Interval<Type>, List<Interval<Type>>> entry : intervals.entrySet()) {
			if(entry.getKey().intersects(target))
				for(Interval<Type> interval : entry.getValue())
					result.add(interval);
			else if(entry.getKey().getStart() > target.getEnd())
				break;
		}

		if(target.getStart() < center && leftNode != null)
			result.addAll(leftNode.query(target));
		if(target.getEnd() > center && rightNode != null)
			result.addAll(rightNode.query(target));
		return result;
	}

	public long getCenter() {
		return center;
	}

	public void setCenter(long center) {
		this.center = center;
	}

	public IntervalNode<Type> getLeft() {
		return leftNode;
	}

	public void setLeft(IntervalNode<Type> left) {
		this.leftNode = left;
	}

	public IntervalNode<Type> getRight() {
		return rightNode;
	}

	public void setRight(IntervalNode<Type> right) {
		this.rightNode = right;
	}

	/**
	 * @param set the set to look on
	 * @return	  the median of the set, not interpolated
	 */
	private Long getMedian(SortedSet<Long> set) {
		int i = 0;
		int middle = set.size() / 2;
		for(Long point : set) {
			if(i == middle)
				return point;
			i++;
		}
		return null;
	}

	@Override
	public String toString() {
		StringBuffer sb = new StringBuffer();
		sb.append(center + ": ");
		for(Entry<Interval<Type>, List<Interval<Type>>> entry : intervals.entrySet()) {
			sb.append("[" + entry.getKey().getStart() + "," + entry.getKey().getEnd() + "]:{");
			for(Interval<Type> interval : entry.getValue()) {
				sb.append("("+interval.getStart()+","+interval.getEnd()+","+interval.getData()+")");
			}
			sb.append("} ");
		}
		return sb.toString();
	}

}
