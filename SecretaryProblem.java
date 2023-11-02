package knapsack;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.ThreadLocalRandom;

public class SecretaryProblem {
	public SecretaryProblem(int n, double w) {
		super();
		this.n = n;
		W = w;
		m= (int) (n / Math.exp(1)); // Number of steps for sampling
	}

	private static  int n;
	private static  double  W;
	private static  int m;
	public static void optimisticSolutionOfSecretartProblem(List<Item> completeItemList) {
		int counter = 0;
		List<Item> samplingItemList = new ArrayList<>();
		List<Item> finalSetOfOptimisticSolution = new ArrayList<>();

		for (Item item : completeItemList) {
			if (counter++ < m) { // Sampling n/e items in order of arrival
				if (samplingItemList.size() < W)
					samplingItemList.add(item);
				else {
					Collections.sort(samplingItemList);
					Collections.reverse(samplingItemList);
					if (item.getUtility() > samplingItemList.get(samplingItemList.size() - 1).getUtility()) {
						samplingItemList.remove(samplingItemList.size() - 1);
						samplingItemList.add(item);

					}
				}
				continue;
			}

			if (evaluateCurrentItemAccordingToOptimisticSecretaryProblem(item, samplingItemList,
					getItemListAccumulativeWeight(finalSetOfOptimisticSolution)) == true) { // returns true if new item
																							// was evaluated good enough
																							// to add to set-solution
				finalSetOfOptimisticSolution.add(item);

			}
		}
		System.out.println("Optimistic solution: " + getItemListAccumulativeUtility(finalSetOfOptimisticSolution));
	}

	public static boolean evaluateCurrentItemAccordingToOptimisticSecretaryProblem(Item contenderItem,
			List<Item> sampledItems, double w) {
		if (contenderItem.getUtility() > sampledItems.get(sampledItems.size() - 1).getUtility()
				&& (w + contenderItem.getWeight()) <= W) { // comparing new item and least valuable item of samplelist

			sampledItems.remove(sampledItems.size() - 1); // In contrast to virtual solution the new item is not added
															// to sampled list

			Collections.sort(sampledItems);
			Collections.reverse(sampledItems);// Is required since sort returns increasing list

			return true;

		}
		return false;
	}

	public void virtualSolutionOfScretaryProblem(List<Item> items, List<Item> samplingItemList,
			List<Item> completeItemList) {
		int counter = 0;

		for(Item item :completeItemList) {	

			if (counter++< m) { // Sampling n/e items in order of arrival
				samplingItemList.add(item);
				if (samplingItemList.size() < W)
					samplingItemList.add(item);
				else {
					Collections.sort(samplingItemList);
					Collections.reverse(samplingItemList);
					if (item.getUtility() > samplingItemList.get(samplingItemList.size() - 1).getUtility()) {
						samplingItemList.remove(samplingItemList.size() - 1);
						samplingItemList.add(item);

					}
				}
				continue;
			}
			// Item added if deemed valuable enough
			if (evaluateCurrentItemAccordingToVirtualSecretaryProblem(item, samplingItemList,
					getItemListAccumulativeWeight(items)) == true) {
				items.add(item);
			}

			if (counter == n) // Stop running when n items reached
				break;
		}

		System.out.println("Virtual algorithm result is:" + getItemListAccumulativeUtility(items));

	}

	public static int calculateOptimalSelection(List<Item> completeItemList) {
		Collections.sort(completeItemList);
		Collections.reverse(completeItemList);
		List<Item> optimalList = new ArrayList<>();

		for (Item item : completeItemList) {
			if (getItemListAccumulativeWeight(optimalList) + item.getWeight() > W)
				break;
			optimalList.add(item);
		}

		return getItemListAccumulativeUtility(optimalList);
	}

	private static boolean evaluateCurrentItemAccordingToVirtualSecretaryProblem(Item contenderItem,
			List<Item> sampledItems, double w) {

		if (contenderItem.getUtility() > sampledItems.get(sampledItems.size() - 1).getUtility()
				&& (w + contenderItem.getWeight()) <= W) { // comparing new item and least valuable item of samplelist

			sampledItems.remove(sampledItems.size() - 1);
			sampledItems.add(contenderItem);

			Collections.sort(sampledItems);
			Collections.reverse(sampledItems);

			return true;

		}

		return false;
	}

	private static Item generateItem(int counter) {// creates item with a weight of 1
		int weight = determineWeight();
		int utility = determineUtility();
		Item item = new Item(counter, weight, utility);
		return item;
	}

	private static int determineUtility() {
		int min = 1;
		int max = 10;
		int utility = ThreadLocalRandom.current().nextInt(min, max + 1); // nextInt is exclusive of top value so add +1
		return utility;
	}

	private static int determineWeight() {
		int min = 1;
		int max = 10;
		int weight = ThreadLocalRandom.current().nextInt(min, max + 1); // nextInt is
		// exclusive of top value so add +1
		return 1;

	}

	private static double getItemListAccumulativeWeight(List<Item> items) {
		double sum = 0;
		if (items.size() == 0)
			return 0;
		for (Item item : items) {
			sum += item.getWeight();
		}
		return sum;
	}

	private static int getItemListAccumulativeUtility(List<Item> items) {
		int sum = 0;
		for (Item item : items) {
			sum += item.getUtility();
		}
		return sum;
	}

	
}


